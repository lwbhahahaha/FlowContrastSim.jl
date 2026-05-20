# ── Contrast agent transport ──
# Plug-flow model: contrast arrives at each segment after cumulative transit time
# with dispersion broadening as it propagates deeper

"""
    ContrastResult

Holds time-resolved contrast concentration.

* When `segment_ids` is empty (legacy behavior), `concentration[s, ti]` is the
  concentration in segment `s` (`1:nseg` of the source tree) at time `ti`.
* When `segment_ids` is non-empty (sparse mode, used for very large trees),
  `concentration[i, ti]` is the concentration in segment `segment_ids[i]`.
  Segments not listed are implicitly zero.

The sparse variant lets us run contrast simulation on 25M+ segment 8 μm trees
without allocating a 40 GB matrix per tree.

`arrival_s` is the physical bolus arrival time in seconds, computed from the
hemodynamics-derived per-segment transit times via topological BFS. Length is
`length(tree.segment_start)` (full tree, NOT filtered) so downstream consumers
can index by segment_id directly. Unreachable / zero-flow segments are `Inf`.
"""
struct ContrastResult
    times::Vector{Float64}                # seconds
    concentration::Matrix{Float64}        # [n_rows x n_timesteps] mg/mL
    outlet_concentration::Matrix{Float64} # [n_rows x n_timesteps] mg/mL
    segment_ids::Vector{Int}              # empty ⇒ dense (row i == segment i)
    arrival_s::Vector{Float64}            # length n_seg full tree; Inf = unreachable
end
ContrastResult(times, C, Cout) = ContrastResult(times, C, Cout, Int[], Float64[])
ContrastResult(times, C, Cout, sids) = ContrastResult(times, C, Cout, sids, Float64[])

# ── Gamma-variate input curve (classic bolus) ──
function gamma_variate_input(times::Vector{Float64};
                             amplitude::Float64=5.0,
                             t0::Float64=0.5,
                             tmax::Float64=4.0,
                             alpha::Float64=3.0)
    C = zeros(length(times))
    for (i, t) in enumerate(times)
        t <= t0 && continue
        tp = (t - t0) / (tmax - t0)
        tp <= 0 && continue
        C[i] = amplitude * tp^alpha * exp(alpha * (1.0 - tp))
    end
    return C
end

# ── Interpolate concentration on uniform time grid ──
function _interp_uniform(C::Vector{Float64}, times::Vector{Float64}, t_query::Float64)
    t_query <= times[1] && return 0.0
    t_query >= times[end] && return C[end]
    dt = times[2] - times[1]
    idx_f = (t_query - times[1]) / dt + 1.0
    i_lo = floor(Int, idx_f)
    i_lo = clamp(i_lo, 1, length(C) - 1)
    frac = idx_f - i_lo
    return C[i_lo] * (1.0 - frac) + C[i_lo + 1] * frac
end

# ── Compute cumulative arrival time from root for each segment ──
function _compute_arrival_times(tree::FlowTree, hemo::HemodynamicsResult)
    nseg = length(tree.segment_start)
    arrival = fill(Inf, nseg)

    order = _topo_order(tree)
    root_segs = Set{Int}()
    for c in tree.children[tree.root_vertex]
        seg = tree.incoming_segment[c]
        seg != 0 && push!(root_segs, seg)
    end

    for s in order
        tau_s = hemo.transit_time_s[s]
        if s in root_segs
            arrival[s] = isfinite(tau_s) ? tau_s / 2.0 : Inf
        else
            start_v = tree.segment_start[s]
            pseg = tree.incoming_segment[start_v]
            if pseg > 0 && isfinite(arrival[pseg])
                parent_tau = hemo.transit_time_s[pseg]
                parent_exit = arrival[pseg] + (isfinite(parent_tau) ? parent_tau / 2.0 : 0.0)
                arrival[s] = parent_exit + (isfinite(tau_s) ? tau_s / 2.0 : 0.0)
            end
        end
    end
    return arrival
end

# ── Plug-flow contrast simulation ──
#
# Physics:
#   - Per-segment arrival[s] is the cumulative transit time along the path from
#     the root to segment s's midpoint, computed from hemodynamics
#     (vol[s] / |flow[s]|) via topological BFS. No log-compression, no
#     calibration — what falls out of Poiseuille + Pries viscosity + the
#     downstream capillary-bed R is what we use.
#   - At segment s, C_s(t) = C_root(t - arrival[s], dispersed). Dispersion
#     stretches the bolus in time by disp_factor and rescales amplitude by
#     1/disp_factor to preserve mass (Gaussian-ish, but applied to the gamma
#     shape — a common in-vivo empirical approximation):
#         disp_factor = sqrt(1 + arrival / t_dispersion_s)
#     The dispersion time scale `t_dispersion_s` is empirical for branching
#     pulsatile vasculature (≈ 1–3 s in vivo coronary; configurable per call /
#     per config).
#   - Unreachable / zero-flow segments keep arrival=Inf ⇒ identically zero
#     contrast (physical: no flow means no contrast delivered).
#
# `max_arrival_s` is accepted as a kwarg for API stability but is IGNORED —
# previously it drove a log-compression of arrival times into a fixed window,
# which was a visualization hack that contaminated the peak-time estimate and
# downstream voxelization. Real physics now governs arrival end-to-end.
function simulate_contrast(tree::FlowTree, hemo::HemodynamicsResult;
                           dt::Float64=0.05,
                           t_end::Float64=30.0,
                           root_input::Union{Nothing, Vector{Float64}}=nothing,
                           amplitude::Float64=5.0,
                           t0::Float64=0.5,
                           tmax::Float64=4.0,
                           alpha::Float64=3.0,
                           max_arrival_s::Float64=0.0,        # accepted but ignored (see note above)
                           t_dispersion_s::Float64=3.0,
                           min_diameter_um::Float64=0.0)
    nseg = length(tree.segment_start)
    times = collect(0.0:dt:t_end)
    nt = length(times)

    # Sparse mode: only simulate segments with diameter ≥ min_diameter_um.
    # Essential for 8 μm capillary trees (25M segs × 200 frames × 8 B = 40 GB
    # dense, but only the visible ≥50 μm skeleton needs time-resolved data).
    segment_ids = Int[]
    if min_diameter_um > 0.0
        thresh_cm = min_diameter_um * 1e-4
        segment_ids = findall(d -> d >= thresh_cm, tree.segment_diameter_cm)
        isempty(segment_ids) && error("No segments with diameter ≥ $(min_diameter_um) μm in $(tree.name)")
    end
    nrows = isempty(segment_ids) ? nseg : length(segment_ids)

    # Physical arrival times (no log-compression, no calibration)
    arrival = _compute_arrival_times(tree, hemo)

    # Gamma-variate input (informational; the per-segment loop samples it directly below)
    C_root = if root_input !== nothing
        root_input
    else
        gamma_variate_input(times; amplitude=amplitude, t0=t0, tmax=tmax, alpha=alpha)
    end

    # Concentration arrays (sized by nrows: either full nseg, or filtered subset)
    C = zeros(nrows, nt)
    C_out = zeros(nrows, nt)

    # Iterate over rows of the output matrix, mapping row → tree segment index.
    @inbounds for i in 1:nrows
        s = isempty(segment_ids) ? i : segment_ids[i]
        a = arrival[s]
        !isfinite(a) && continue

        disp_factor = sqrt(1.0 + a / t_dispersion_s)
        for ti in 1:nt
            t = times[ti]
            t_shifted = t - a
            t_shifted <= 0 && continue

            # Sample dispersed bolus: time axis stretched by disp_factor,
            # amplitude rescaled by 1/disp_factor (preserves area under curve).
            t_input = t0 + (t_shifted - t0) / disp_factor
            if t_input > t0
                tp = (t_input - t0) / (tmax - t0)
                if tp > 0
                    c_val = amplitude * tp^alpha * exp(alpha * (1.0 - tp)) / disp_factor
                    C[i, ti] = max(c_val, 0.0)
                    C_out[i, ti] = C[i, ti]
                end
            end
        end
    end

    return ContrastResult(times, C, C_out, segment_ids, arrival)
end

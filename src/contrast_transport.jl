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
"""
struct ContrastResult
    times::Vector{Float64}                # seconds
    concentration::Matrix{Float64}        # [n_rows x n_timesteps] mg/mL
    outlet_concentration::Matrix{Float64} # [n_rows x n_timesteps] mg/mL
    segment_ids::Vector{Int}              # empty ⇒ dense (row i == segment i)
end
ContrastResult(times, C, Cout) = ContrastResult(times, C, Cout, Int[])

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

# ── Plug-flow contrast simulation with time compression for visualization ──
function simulate_contrast(tree::FlowTree, hemo::HemodynamicsResult;
                           dt::Float64=0.05,
                           t_end::Float64=30.0,
                           root_input::Union{Nothing, Vector{Float64}}=nothing,
                           amplitude::Float64=5.0,
                           t0::Float64=0.5,
                           tmax::Float64=4.0,
                           alpha::Float64=3.0,
                           max_arrival_s::Float64=0.0,
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

    # Compute raw arrival times
    arrival_raw = _compute_arrival_times(tree, hemo)

    # Compress arrival times: map [0, p99_arrival] -> [0, effective_window]
    finite_arr = sort(arrival_raw[isfinite.(arrival_raw)])
    if isempty(finite_arr)
        return ContrastResult(times, zeros(nrows, nt), zeros(nrows, nt), segment_ids)
    end

    # Auto-determine max arrival if not specified
    if max_arrival_s <= 0
        # Use t_end - bolus_peak as the window for contrast to fill tree
        max_arrival_s = t_end - tmax
    end

    # Use p99 as the reference for compression
    p99_raw = length(finite_arr) > 10 ? finite_arr[min(length(finite_arr), Int(ceil(0.99 * length(finite_arr))))] : finite_arr[end]
    p99_raw = max(p99_raw, 1.0)

    # Log-compress: arrival_compressed = max_arrival * log(1 + arrival_raw/scale) / log(1 + p99/scale)
    # scale controls how much compression (smaller = more compression of large values)
    scale = max(finite_arr[max(1, length(finite_arr) / 4 |> x -> Int(ceil(x)))], 0.1)  # p25 as scale
    log_norm = log(1.0 + p99_raw / scale)

    arrival = fill(Inf, nseg)
    for s in 1:nseg
        isfinite(arrival_raw[s]) || continue
        arrival[s] = max_arrival_s * log(1.0 + arrival_raw[s] / scale) / log_norm
    end

    # Gamma-variate input
    C_root = if root_input !== nothing
        root_input
    else
        gamma_variate_input(times; amplitude=amplitude, t0=t0, tmax=tmax, alpha=alpha)
    end

    # Dispersion broadening
    t_disp = 3.0  # seconds

    # Concentration arrays (sized by nrows: either full nseg, or filtered subset)
    C = zeros(nrows, nt)
    C_out = zeros(nrows, nt)

    # Iterate over rows of the output matrix, mapping row → tree segment index.
    @inbounds for i in 1:nrows
        s = isempty(segment_ids) ? i : segment_ids[i]
        !isfinite(arrival[s]) && continue

        for ti in 1:nt
            t = times[ti]
            t_shifted = t - arrival[s]
            t_shifted <= 0 && continue

            # Apply dispersion: broader bolus for deeper segments
            disp_factor = sqrt(1.0 + arrival[s] / t_disp)

            # Sample dispersed input curve
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

    return ContrastResult(times, C, C_out, segment_ids)
end

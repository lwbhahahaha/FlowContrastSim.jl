# ── Contrast agent transport on the tree ──
#
# Two solution modes, selected by whether an `aif::Vector{Float64}` is passed:
#
# 1. AIF-driven (Taylor-Aris analytical PDE solution). The 1D advection-
#    dispersion equation
#        ∂C/∂t + v ∂C/∂x = D_eff ∂²C/∂x²
#    is linear with steady flow + linear branching. Its exact path-integrated
#    Green's function solution at segment s, sampled at midpoint, is
#        C_s(t) = (AIF ⊛ G_{μ_s, σ_s})(t)
#    with
#        μ_s    = Σ_{ancestors} τ_a   + τ_s / 2          (cumulative mean transit time)
#        σ_s²   = Σ_{ancestors} σ_a²  + σ_seg_s² / 2    (cumulative Taylor variance)
#        σ_seg² = 2 D_eff(R,v) × L / v³
#        D_eff  = D_mol + (R v)² / (48 D_mol)             (Taylor-Aris)
#    All inputs come from anatomy + hemodynamics + the iodine molecular
#    diffusivity D_mol (≈ 1.5e-9 m²/s for iodinated contrast in plasma).
#    No free knobs.
#
# 2. Legacy gamma-variate (kept for backward compatibility with existing
#    configs and the pre-protocol test suite). Uses an empirical
#        disp_factor = sqrt(1 + arrival_s / t_dispersion_s)
#    with `t_dispersion_s` a configurable scalar. Not derived from segment
#    geometry; do not use for new physics work.
#
# Mode 1 activates when the caller passes a non-`nothing` `aif`. Both modes
# share the same `arrival_s` BFS skeleton, the sparse diameter-filter mask,
# and the `ContrastResult` struct shape.

"""
    ContrastResult

Holds time-resolved contrast concentration.

* When `segment_ids` is empty (legacy behavior), `concentration[s, ti]` is the
  concentration in segment `s` (`1:nseg` of the source tree) at time `ti`.
* When `segment_ids` is non-empty (sparse mode, used for very large trees),
  `concentration[i, ti]` is the concentration in segment `segment_ids[i]`.
  Segments not listed are implicitly zero.

The sparse variant lets us run contrast simulation on 25M+ segment 8 μm trees
without allocating a 40 GB matrix per tree. When the `peak_time_s` kwarg of
`simulate_contrast` is set, `concentration` has a single column (one snapshot
at the requested acquisition time) and `times` has length 1.

`arrival_s` is the physical bolus arrival time in seconds, computed from the
hemodynamics-derived per-segment transit times via topological BFS. Length is
`length(tree.segment_start)` (full tree, NOT filtered) so downstream consumers
can index by segment_id directly. Unreachable / zero-flow segments are `Inf`.

`arrival_variance_s2` is the per-segment cumulative Taylor-Aris variance σ²
of bolus arrival time at the segment midpoint, in s². Populated only on the
PDE / AIF path; on the legacy gamma-variate path it is `Inf` for every
segment. Length is `length(tree.segment_start)`.
"""
struct ContrastResult
    times::Vector{Float64}                # seconds
    concentration::Matrix{Float64}        # [n_rows x n_timesteps] mg/mL
    outlet_concentration::Matrix{Float64} # [n_rows x n_timesteps] mg/mL
    segment_ids::Vector{Int}              # empty ⇒ dense (row i == segment i)
    arrival_s::Vector{Float64}            # length n_seg full tree; Inf = unreachable
    arrival_variance_s2::Vector{Float64}  # length n_seg full tree; Inf when PDE path inactive
end
ContrastResult(times, C, Cout) = ContrastResult(times, C, Cout, Int[], Float64[], Float64[])
ContrastResult(times, C, Cout, sids) = ContrastResult(times, C, Cout, sids, Float64[], Float64[])
ContrastResult(times, C, Cout, sids, arr) = ContrastResult(times, C, Cout, sids, arr, fill(Inf, length(arr)))

# Molecular diffusivity of iodinated contrast in plasma (m²/s).
# Iohexol / iomeprol literature: ~1.3–1.6e-9 at 37 °C. We use 1.5e-9 as a
# round midpoint. Sensitivity is weak in the advection-dominated regime
# (Pe ≫ 1 everywhere ≥ capillary) because the Taylor term R²v²/(48 D_mol)
# scales as 1/D_mol while the molecular term D_mol vanishes — net dependence
# on D_mol is logarithmic, not first-order.
const D_MOL_PLASMA_DEFAULT_M2_S = 1.5e-9

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

# ── Taylor-Aris dispersion: per-segment time-domain variance σ² ──
#
# For 1D Poiseuille flow in a tube, the effective axial dispersion coefficient is
#     D_eff = D_mol + (R v)² / (48 D_mol)        (Taylor 1953, Aris 1956)
# The bolus arrival time at the outlet of a single segment of length L follows
# a Gaussian with mean L/v and variance
#     σ²_t = 2 D_eff L / v³.
# Variance is additive along independent serial segments (variance of a sum
# of independent transit times = sum of variances), giving the cumulative
# variance up to any downstream segment as a pure BFS accumulation — the same
# topology pass used for mean arrival time.
"""
    _segment_taylor_variance(R_m, v_m_s, L_m; D_mol)

Per-segment Taylor-Aris time-variance σ² (s²) for a bolus traversing a tube
from inlet to outlet. All arguments in SI units. Returns 0 for degenerate
geometry (zero flow / zero length / zero radius) — that segment contributes
no dispersion to any downstream path.
"""
function _segment_taylor_variance(R_m::Float64, v_m_s::Float64, L_m::Float64;
                                   D_mol::Float64=D_MOL_PLASMA_DEFAULT_M2_S)
    (R_m <= 0.0 || v_m_s <= 0.0 || L_m <= 0.0) && return 0.0
    D_eff = D_mol + (R_m * v_m_s)^2 / (48.0 * D_mol)
    return 2.0 * D_eff * L_m / v_m_s^3
end

# ── Per-segment (arrival mean, arrival variance) at midpoint ──
#
# Topological BFS propagates both moments together. For root segments
# (no parent), the midpoint accumulates half of that segment's transit
# time and half of its Taylor variance. For non-root segments, we ADD the
# parent's full segment variance (inlet→outlet of parent) plus half of
# our own (we sample at our midpoint).
function _compute_arrival_moments(tree::FlowTree, hemo::HemodynamicsResult;
                                   D_mol::Float64=D_MOL_PLASMA_DEFAULT_M2_S)
    nseg = length(tree.segment_start)
    # Pre-compute per-segment inlet→outlet variance once.
    seg_full_var = zeros(nseg)
    @inbounds for s in 1:nseg
        a = tree.vertices[tree.segment_start[s]]
        b = tree.vertices[tree.segment_end[s]]
        L_m   = norm(b - a) * 0.01            # cm → m
        R_m   = tree.segment_diameter_cm[s] * 0.5 * 0.01
        Q_m3s = abs(hemo.segment_flow[s])
        if Q_m3s > 1e-30 && R_m > 0.0 && L_m > 0.0
            v_m_s = Q_m3s / (pi * R_m^2)
            seg_full_var[s] = _segment_taylor_variance(R_m, v_m_s, L_m; D_mol=D_mol)
        end
    end

    arrival_mean = fill(Inf, nseg)
    arrival_var  = fill(Inf, nseg)
    order = _topo_order(tree)
    root_segs = Set{Int}()
    for c in tree.children[tree.root_vertex]
        seg = tree.incoming_segment[c]
        seg != 0 && push!(root_segs, seg)
    end

    for s in order
        tau_s = hemo.transit_time_s[s]
        (!isfinite(tau_s) || tau_s <= 0.0) && continue
        sv = seg_full_var[s]

        if s in root_segs
            arrival_mean[s] = tau_s / 2.0
            arrival_var[s]  = sv / 2.0
        else
            start_v = tree.segment_start[s]
            pseg = tree.incoming_segment[start_v]
            if pseg > 0 && isfinite(arrival_mean[pseg])
                ptau = hemo.transit_time_s[pseg]
                pvar = seg_full_var[pseg]
                arrival_mean[s] = arrival_mean[pseg] +
                                  (isfinite(ptau) ? ptau / 2.0 : 0.0) +
                                  tau_s / 2.0
                arrival_var[s]  = arrival_var[pseg] + pvar / 2.0 + sv / 2.0
            end
        end
    end
    return arrival_mean, arrival_var
end

# ── Discrete Gaussian convolution of AIF with width σ on the same time grid ──
#
# Returns F such that F[i] ≈ (AIF ⊛ G_σ)(times[i]). When σ ≤ dt/2 the
# Gaussian narrower than the time discretization, so we return AIF unchanged
# (delta-function limit; discretization noise would otherwise dominate).
# The kernel is normalized so that its discrete sum × dt = 1 — preserves
# the integral of AIF exactly (modulo window truncation at ±4σ, < 0.01 %).
function _gaussian_convolve(aif::Vector{Float64}, dt::Float64, σ::Float64)
    nt = length(aif)
    if σ <= 0.5 * dt
        return copy(aif)
    end
    hw = max(1, ceil(Int, 4.0 * σ / dt))
    kernel = Vector{Float64}(undef, 2 * hw + 1)
    @inbounds for k in -hw:hw
        kernel[k + hw + 1] = exp(-(k * dt)^2 / (2.0 * σ * σ))
    end
    Z = sum(kernel) * dt
    Z > 0.0 && (kernel ./= Z)

    F = zeros(nt)
    @inbounds for i in 1:nt
        acc = 0.0
        for m in -hw:hw
            j = i + m
            if 1 <= j <= nt
                acc += aif[j] * kernel[m + hw + 1]
            end
        end
        F[i] = acc * dt
    end
    return F
end

# Linear-interpolation sample of F (already on `times` grid) at time `t`.
# Returns 0 outside [0, t_end]; F is the convolved AIF so it should be
# essentially zero outside the AIF support window.
function _sample_grid(F::Vector{Float64}, dt::Float64, t::Float64)
    t < 0.0 && return 0.0
    nt = length(F)
    idx_f = t / dt + 1.0
    i_lo = floor(Int, idx_f)
    i_lo < 1 && return 0.0
    i_lo >= nt && return 0.0
    frac = idx_f - i_lo
    return F[i_lo] * (1.0 - frac) + F[i_lo + 1] * frac
end

# ── Single-point Gaussian convolution: (AIF ⊛ G_σ)(t) for one t value ──
#
# Equivalent to `_gaussian_convolve` followed by `_sample_grid`, but
# computes only the one query point. Used for peak-only mode: per
# segment we need C(s, t*) at one acquisition time, not the whole curve.
#
# Discrete convolution with normalization: the window is renormalized so
# that the discrete kernel sum × dt = 1 within the window. This handles
# truncation at the AIF support boundary (when t is near 0 or t_end) as
# well as inside-window evaluations consistently.
function _gaussian_convolve_at(aif::Vector{Float64}, dt::Float64, σ::Float64, t::Float64)
    t < 0.0 && return 0.0
    nt = length(aif)
    if σ <= 0.5 * dt
        # Delta limit: just linearly interpolate AIF at t
        idx_f = t / dt + 1.0
        i_lo = floor(Int, idx_f)
        i_lo < 1 && return 0.0
        i_lo >= nt && return 0.0
        frac = idx_f - i_lo
        return aif[i_lo] * (1.0 - frac) + aif[i_lo + 1] * frac
    end

    hw = max(1, ceil(Int, 4.0 * σ / dt))
    i_center = round(Int, t / dt) + 1
    j_lo = max(1, i_center - hw)
    j_hi = min(nt, i_center + hw)
    inv_2σ² = 1.0 / (2.0 * σ * σ)

    sum_aw = 0.0
    sum_w  = 0.0
    @inbounds for j in j_lo:j_hi
        t_off = (j - 1) * dt - t
        w = exp(-t_off * t_off * inv_2σ²)
        sum_aw += aif[j] * w
        sum_w  += w
    end
    return sum_w > 0.0 ? sum_aw / sum_w : 0.0
end

# ── Contrast transport: AIF (Taylor-Aris PDE) or legacy gamma-variate ──
"""
    simulate_contrast(tree, hemo; dt, t_end, aif=nothing, D_mol_m2_s=1.5e-9,
                       peak_time_s=nothing, min_diameter_um=0, ...) -> ContrastResult

Per-segment contrast concentration time series (or single snapshot).

When `aif` is provided (a vector aligned with `0:dt:t_end`), the per-segment
concentration is the analytical Green's function solution of the 1D
advection-dispersion equation on the tree: for each segment s,

    C_s(t) = (AIF ⊛ G_{μ_s, σ_s})(t)

with `μ_s, σ_s` accumulated from per-segment geometry + flow via Taylor-Aris.
`D_mol_m2_s` (default 1.5e-9) is the iodine molecular diffusivity.

When `aif === nothing`, the legacy hand-tuned gamma-variate input is used
with the empirical `disp_factor = sqrt(1 + arrival/t_dispersion_s)`.

Sparse-mode kwarg `min_diameter_um > 0` restricts simulation to segments
above that diameter (saves O(N_seg × N_t) memory on 8 μm capillary trees).

When `peak_time_s::Float64` is set (AIF path only), only the single
snapshot `C(s, peak_time_s)` is computed and stored. The returned
`ContrastResult.concentration` is `nrows × 1` and `.times = [peak_time_s]`.
This is the mode used by SVP-style perfusion (one V2 acquisition at a
clinician-chosen post-injection time); it cuts memory and compute by
~N_t× compared to the full time series.
"""
function simulate_contrast(tree::FlowTree, hemo::HemodynamicsResult;
                           dt::Float64=0.05,
                           t_end::Float64=30.0,
                           aif::Union{Nothing, Vector{Float64}}=nothing,
                           D_mol_m2_s::Float64=D_MOL_PLASMA_DEFAULT_M2_S,
                           peak_time_s::Union{Nothing, Float64}=nothing,
                           # Legacy gamma-variate parameters; used only when `aif === nothing`.
                           amplitude::Float64=5.0,
                           t0::Float64=0.5,
                           tmax::Float64=4.0,
                           alpha::Float64=3.0,
                           max_arrival_s::Float64=0.0,        # accepted but ignored
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

    using_aif = aif !== nothing
    if using_aif
        length(aif) == nt || error("aif length $(length(aif)) ≠ times length $(nt); regenerate AIF with the same dt and t_end")

        # ── Taylor-Aris PDE path ──
        arrival_mean, arrival_var = _compute_arrival_moments(tree, hemo; D_mol=D_mol_m2_s)

        # Peak-only snapshot mode: just one time point per segment.
        if peak_time_s !== nothing
            t_peak = Float64(peak_time_s)
            (t_peak < 0.0 || t_peak > times[end]) &&
                error("peak_time_s=$(t_peak) outside simulation window [0, $(times[end])]")
            out_times = [t_peak]
            C = zeros(nrows, 1)
            @inbounds for i in 1:nrows
                s = isempty(segment_ids) ? i : segment_ids[i]
                μ = arrival_mean[s]
                v = arrival_var[s]
                (!isfinite(μ) || !isfinite(v)) && continue
                σ = sqrt(max(v, 0.0))
                c_val = _gaussian_convolve_at(aif, dt, σ, t_peak - μ)
                if c_val > 0.0
                    C[i, 1] = c_val
                end
            end
            return ContrastResult(out_times, C, copy(C), segment_ids, arrival_mean, arrival_var)
        end

        # Full time-series mode.
        C = zeros(nrows, nt)
        C_out = zeros(nrows, nt)
        @inbounds for i in 1:nrows
            s = isempty(segment_ids) ? i : segment_ids[i]
            μ = arrival_mean[s]
            v = arrival_var[s]
            (!isfinite(μ) || !isfinite(v)) && continue

            σ = sqrt(max(v, 0.0))
            F = _gaussian_convolve(aif, dt, σ)

            for ti in 1:nt
                t_local = times[ti] - μ
                c_val = _sample_grid(F, dt, t_local)
                if c_val > 0.0
                    C[i, ti] = c_val
                    C_out[i, ti] = c_val
                end
            end
        end

        return ContrastResult(times, C, C_out, segment_ids, arrival_mean, arrival_var)
    end

    # ── Legacy gamma-variate path (empirical sqrt dispersion) ──
    peak_time_s === nothing || error("peak_time_s is only supported on the AIF / PDE path; supply `aif=...`")
    arrival = _compute_arrival_times(tree, hemo)
    C = zeros(nrows, nt)
    C_out = zeros(nrows, nt)

    @inbounds for i in 1:nrows
        s = isempty(segment_ids) ? i : segment_ids[i]
        a = arrival[s]
        !isfinite(a) && continue

        disp_factor = sqrt(1.0 + a / t_dispersion_s)
        for ti in 1:nt
            t = times[ti]
            t_shifted = t - a
            t_shifted <= 0 && continue

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

#!/usr/bin/env julia
# LAD proximal-stenosis sweep with closed-loop autoregulation (analytical).
#
# Models the at-rest state with active autoregulatory tone in the resistance
# vessel band [8, 400] μm. As proximal stenosis is applied, the in-band
# arterioles dilate uniformly to compensate, up to a 1.6× reserve cap
# (Wong-Molloi 2008 empirical factor). When the reserve is exhausted, root
# flow falls below the no-stenosis target.
#
# Tree: max-dilated CSV. Apply per-segment factor f ∈ [F_MIN, F_MAX] to in-
# band segments. f = F_MIN = 0.625 → at-rest geometry, f = F_MAX = 1.0 →
# max-dilated geometry.
#
# Closed-form f estimator (avoids ~80 hemo bisection calls):
#
#   R_total(s, f)  =  R_other(s)  +  b / f^4
#
#       b           = Poiseuille R contribution of the arteriole band,
#                     evaluated at f=1.0 (max-dilated). Independent of s
#                     because stenosis is on dias_lad1 (>400 μm, not in band).
#       R_other(s)  = everything else (capillary + cap-bed + conduit + the
#                     stenosed dias_lad1 segments). Depends on s.
#
#   ΔP / Q(s, f) = R_other(s) + b / f^4
#
# Pre-compute b once from two no-stenosis hemos (Q at F_MIN and Q at F_MAX,
# s=0). Then per stenosis level, one Q(s, F_MAX) measurement gives R_other(s),
# and the required f solves Q(s, f) = Q_target:
#
#   f^4 = b / (ΔP/Q_target − R_other(s))
#
# One verification hemo per s; one Newton refinement if |Q − Q_target|/Q_target
# > 2 %. Total: 2 + ~30 hemos for the sweep ≈ 40 min instead of 90.
#
# The model ignores Pries viscosity's mild dependence on f (small arterioles
# change diameter, slightly shifting their η_rel). Newton refinement absorbs
# that residual.
#
# Usage:
#   julia --project=. scripts/lad_stenosis_autoreg_sweep.jl \
#         <lad_csv_max_dilated> <baseline_config.toml> [output_csv]

using FlowContrastSim
using Printf
using Statistics
import FlowContrastSim: load_tree, compute_hemodynamics, load_flow_config

length(ARGS) >= 2 || error("Usage: lad_stenosis_autoreg_sweep.jl <lad_csv> <baseline_config> [output_csv]")

const CSV_PATH    = ARGS[1]
const CONFIG_PATH = ARGS[2]
const cfg_base    = splitext(basename(CONFIG_PATH))[1]
const OUT_CSV     = length(ARGS) >= 3 ? ARGS[3] : "scripts/lad_stenosis_autoreg_$(cfg_base).csv"

const PROXIMAL_LABEL = "dias_lad1"
const STENOSIS_PCTS  = 0.0:10.0:90.0

const D_LOW_UM       = 8.0       # Wong-Molloi resistance band lower
const D_HIGH_UM      = 400.0     # Wong-Molloi resistance band upper
const F_MIN          = 0.625     # at-rest (full autoregulatory tone)
const F_MAX          = 1.0       # max dilated (1.6× reserve fully used)
const REFINE_TOL     = 0.02      # Newton refinement when |Q-Q_target|/Q_target > this

cfg = load_flow_config(CONFIG_PATH)
mass_lad  = get(cfg.territory_masses_g, "LAD", 0.0)
ROOT_P_MMHG = cfg.root_pressure_mmhg
TERM_P_MMHG = cfg.terminal_pressure_mmhg
ROOT_P_PA   = ROOT_P_MMHG * 133.322
TERM_P_PA   = TERM_P_MMHG * 133.322
DP_MMHG     = ROOT_P_MMHG - TERM_P_MMHG

println("LAD stenosis sweep WITH autoregulation (analytical solver)")
println("  tree CSV: $(CSV_PATH)")
println("  config:   $(CONFIG_PATH)")
println("  output:   $(OUT_CSV)")
println("  proximal stenosis label: $(PROXIMAL_LABEL)")
println("  autoreg band: [$D_LOW_UM, $D_HIGH_UM] μm")
println("  dilation factor: f ∈ [$F_MIN, $F_MAX]   (reserve = $(round(1/F_MIN; digits=2))x)")
println("  BCs: $(ROOT_P_MMHG) → $(TERM_P_MMHG) mmHg, Hct=$(cfg.discharge_hematocrit)")
println("  cap_R=$(cfg.capillary_bed_R_per_100g_mmHgmin_ml) mmHg·min/mL/100g, LAD mass=$(mass_lad) g")
println("-" ^ 80)

print("loading LAD tree ... "); flush(stdout)
t0 = time()
tree = load_tree("LAD", CSV_PATH)
@printf("done in %.1fs, %d segs, %d verts\n", time()-t0, length(tree.segment_start), length(tree.vertices))

in_band_idx = findall(i -> D_LOW_UM <= tree.segment_diameter_cm[i] * 1e4 <= D_HIGH_UM,
                     eachindex(tree.segment_diameter_cm))
prox_idx    = findall(==(PROXIMAL_LABEL), tree.segment_label)
@printf("  autoreg-band segs: %d  (%.1f%%)\n",
        length(in_band_idx), 100*length(in_band_idx)/length(tree.segment_diameter_cm))
@printf("  proximal LAD segs: %d   d range [%.0f, %.0f] μm\n",
        length(prox_idx),
        minimum(tree.segment_diameter_cm[i]*1e4 for i in prox_idx),
        maximum(tree.segment_diameter_cm[i]*1e4 for i in prox_idx))
flush(stdout)

orig_d_cm = copy(tree.segment_diameter_cm)

function hemo_root_flow!(tree, orig_d_cm, in_band_idx, prox_idx, f::Float64, s::Float64, cfg, mass_lad)
    @inbounds for i in eachindex(tree.segment_diameter_cm)
        tree.segment_diameter_cm[i] = orig_d_cm[i]
    end
    @inbounds for i in in_band_idx
        tree.segment_diameter_cm[i] = orig_d_cm[i] * f
    end
    scale = 1.0 - s / 100.0
    @inbounds for i in prox_idx
        tree.segment_diameter_cm[i] = orig_d_cm[i] * scale
    end
    hemo = compute_hemodynamics(tree;
        root_pressure=ROOT_P_PA,
        terminal_pressure=TERM_P_PA,
        hematocrit=cfg.discharge_hematocrit,
        capillary_bed_R_per_100g_mmHgmin_ml=cfg.capillary_bed_R_per_100g_mmHgmin_ml,
        territory_mass_g=mass_lad)
    rf = 0.0
    for c in tree.children[tree.root_vertex]
        seg = tree.incoming_segment[c]
        seg != 0 && (rf += hemo.segment_flow[seg])
    end
    return rf * 60e6   # mL/min
end

# ── Step 1: Q_target = Q(s=0, f=F_MIN), and Q(s=0, F_MAX) for model fit ──
print("computing Q_target (s=0, f=F_MIN=$F_MIN) ... "); flush(stdout)
t1 = time()
Q_target = hemo_root_flow!(tree, orig_d_cm, in_band_idx, prox_idx, F_MIN, 0.0, cfg, mass_lad)
@printf("done in %.1fs  →  Q_target = %.2f mL/min\n", time()-t1, Q_target)

print("computing Q(s=0, f=F_MAX=$F_MAX) ... "); flush(stdout)
t1 = time()
Q_s0_FMAX = hemo_root_flow!(tree, orig_d_cm, in_band_idx, prox_idx, F_MAX, 0.0, cfg, mass_lad)
@printf("done in %.1fs  →  Q = %.2f mL/min\n", time()-t1, Q_s0_FMAX)

# Solve for b (band Poiseuille coefficient) using s=0 endpoints:
#   ΔP/Q_target  = R_other(0) + b / F_MIN^4
#   ΔP/Q_s0_FMAX = R_other(0) + b / F_MAX^4   (= R_other(0) + b, since F_MAX=1)
# Subtract:
#   ΔP × (1/Q_target − 1/Q_s0_FMAX) = b × (1/F_MIN^4 − 1)
R_TARGET   = DP_MMHG / Q_target          # mmHg·min/mL
R_S0_FMAX  = DP_MMHG / Q_s0_FMAX
b          = (R_TARGET - R_S0_FMAX) / (1.0/F_MIN^4 - 1.0)
R_other_0  = R_S0_FMAX - b
@printf("model: b = %.4f mmHg·min/mL ; R_other(s=0) = %.4f\n", b, R_other_0)
println("-" ^ 80)

results = Tuple{Float64, Float64, Float64, Float64, Int, Bool}[]
# (stenosis %, prox_d_um, root_flow_mlmin, f_used, n_hemo_calls, reserve_exhausted)

for s in STENOSIS_PCTS
    n_calls = 0
    prox_d_um = mean(orig_d_cm[i] for i in prox_idx) * 1e4 * (1.0 - s/100.0)

    if s == 0.0
        # We already have Q_target.
        push!(results, (0.0, prox_d_um, Q_target, F_MIN, 0, false))
        @printf("stenosis=%5.1f%%  proximal_d=%7.1f μm  root_flow=%8.2f mL/min  f=%.4f  exhausted=false  (n_hemo=0)\n",
                0.0, prox_d_um, Q_target, F_MIN)
        flush(stdout)
        continue
    end

    # 1. Max dilation: gives R_other(s).
    t_h = time()
    Q_s_FMAX = hemo_root_flow!(tree, orig_d_cm, in_band_idx, prox_idx, F_MAX, s, cfg, mass_lad)
    n_calls += 1
    @printf("  s=%.0f%%  Q@F_MAX=%.2f mL/min  (%.1fs)\n", s, Q_s_FMAX, time()-t_h)
    flush(stdout)

    if Q_s_FMAX < Q_target
        # Reserve exhausted: report at f = F_MAX.
        push!(results, (s, prox_d_um, Q_s_FMAX, F_MAX, n_calls, true))
        @printf("stenosis=%5.1f%%  proximal_d=%7.1f μm  root_flow=%8.2f mL/min  f=%.4f  exhausted=true  (n_hemo=%d)\n",
                s, prox_d_um, Q_s_FMAX, F_MAX, n_calls)
        flush(stdout)
        continue
    end

    # 2. Solve for f analytically:  f^4 = b / (R_TARGET − R_other(s))
    R_s_FMAX  = DP_MMHG / Q_s_FMAX
    R_other_s = R_s_FMAX - b
    denom     = R_TARGET - R_other_s
    if denom <= 0
        # Numerical floor — model says no f can match. Treat as exhausted.
        push!(results, (s, prox_d_um, Q_s_FMAX, F_MAX, n_calls, true))
        @printf("stenosis=%5.1f%%  proximal_d=%7.1f μm  root_flow=%8.2f mL/min  f=%.4f  exhausted=true (R model)  (n_hemo=%d)\n",
                s, prox_d_um, Q_s_FMAX, F_MAX, n_calls)
        flush(stdout)
        continue
    end
    f_est = clamp((b / denom)^0.25, F_MIN, F_MAX)
    @printf("  s=%.0f%%  R_other(s)=%.4f  → analytical f=%.4f\n", s, R_other_s, f_est)
    flush(stdout)

    # 3. Verify.
    t_h = time()
    Q_est = hemo_root_flow!(tree, orig_d_cm, in_band_idx, prox_idx, f_est, s, cfg, mass_lad)
    n_calls += 1
    rel_err = (Q_est - Q_target) / Q_target
    @printf("  s=%.0f%%  Q@analytical = %.2f mL/min (rel err %+.2f%%)  (%.1fs)\n",
            s, Q_est, 100*rel_err, time()-t_h)
    flush(stdout)

    # 4. Optional Newton refinement: dQ/df rises super-linearly in f. Use the
    # b/f^4 model to predict the f_new. R_total = b/f^4 + R_other → Q = ΔP /
    # (b/f^4 + R_other). For a refined estimate, infer effective b from
    # current point and re-solve.
    if abs(rel_err) > REFINE_TOL
        R_now   = DP_MMHG / Q_est
        b_eff   = (R_now - R_other_s) * f_est^4         # back out effective b at this f
        denom2  = R_TARGET - R_other_s
        f_new   = clamp((b_eff / denom2)^0.25, F_MIN, F_MAX)
        t_h = time()
        Q_new = hemo_root_flow!(tree, orig_d_cm, in_band_idx, prox_idx, f_new, s, cfg, mass_lad)
        n_calls += 1
        rel_err2 = (Q_new - Q_target) / Q_target
        @printf("  s=%.0f%%  Newton refine f=%.4f → Q=%.2f (rel err %+.2f%%)  (%.1fs)\n",
                s, f_new, Q_new, 100*rel_err2, time()-t_h)
        flush(stdout)
        if abs(rel_err2) < abs(rel_err)
            f_est = f_new
            Q_est = Q_new
        end
    end

    push!(results, (s, prox_d_um, Q_est, f_est, n_calls, false))
    @printf("stenosis=%5.1f%%  proximal_d=%7.1f μm  root_flow=%8.2f mL/min  f=%.4f  exhausted=false  (n_hemo=%d)\n",
            s, prox_d_um, Q_est, f_est, n_calls)
    flush(stdout)
end

@inbounds for i in eachindex(tree.segment_diameter_cm)
    tree.segment_diameter_cm[i] = orig_d_cm[i]
end

mkpath(dirname(OUT_CSV))
open(OUT_CSV, "w") do io
    println(io, "stenosis_pct,proximal_diameter_um,root_flow_mlmin,autoreg_f,reserve_exhausted")
    for (s, d, q, f, _, ex) in results
        @printf(io, "%.1f,%.4f,%.6f,%.4f,%d\n", s, d, q, f, ex ? 1 : 0)
    end
end
println("-" ^ 80)
println("wrote $(OUT_CSV)")

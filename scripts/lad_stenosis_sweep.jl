#!/usr/bin/env julia
# LAD stenosis sweep: shrink the proximal LAD segments (label == "dias_lad1")
# by 0%, 10%, 20%, ..., 90% and report root flow at each step. Emits a CSV
# at `scripts/lad_stenosis.csv` suitable for a stenosis-vs-flow chart.
#
# Design notes:
#   - Stenosis is applied uniformly to every segment labeled "dias_lad1"
#     (the proximal LAD, 23 consecutive ~0.5 mm segments ≈ 12 mm total).
#     These are the segments where hemodynamically significant stenosis
#     is most often clinically observed.
#   - Diameter of each affected segment: d_new = d_orig × (1 − stenosis/100).
#   - Boundary conditions and viscosity model are unchanged (Pries 1992 +
#     Secomb 2005 ESL). Pressure gradient held at 100 → 15 mmHg.
#   - Tree is loaded once; for each stenosis level we restore the original
#     diameter vector, multiply the proximal indices, and re-solve Poiseuille.
#     Avoids re-parsing the 8 GB CSV 10×.
#
# Usage:
#   julia --project=. scripts/lad_stenosis_sweep.jl [tree_csv_path]
# Default CSV: ../VascularTreeSim.jl/output/lad_segments.csv

using FlowContrastSim
using Printf
using Statistics
import FlowContrastSim: load_tree, compute_hemodynamics

const CSV_PATH = length(ARGS) >= 1 ? ARGS[1] :
    "../VascularTreeSim.jl/output/lad_segments.csv"
const OUT_CSV = "scripts/lad_stenosis.csv"

const ROOT_P_PA = 100.0 * 133.322
const TERM_P_PA = 15.0 * 133.322
const HCT = 0.45
const PROXIMAL_LABEL = "dias_lad1"
const STENOSIS_PCTS = 0.0:10.0:90.0

println("LAD stenosis sweep")
println("  tree CSV: $(CSV_PATH)")
println("  proximal label: $(PROXIMAL_LABEL)")
println("  BCs: $(ROOT_P_PA/133.322) → $(TERM_P_PA/133.322) mmHg, Hct=$(HCT)")
println("-" ^ 80)

print("loading LAD tree ... "); flush(stdout)
t0 = time()
tree = load_tree("LAD", CSV_PATH)
@printf("done in %.1fs, %d segs, %d verts\n", time()-t0, length(tree.segment_start), length(tree.vertices))

# Identify proximal indices
prox_idx = findall(==(PROXIMAL_LABEL), tree.segment_label)
if isempty(prox_idx)
    error("No segments with label '$(PROXIMAL_LABEL)' found in tree.")
end
orig_prox_d_um = [tree.segment_diameter_cm[i] * 1e4 for i in prox_idx]
@printf("  found %d proximal segments, diameter range [%.1f, %.1f] μm\n",
    length(prox_idx),
    minimum(orig_prox_d_um), maximum(orig_prox_d_um))

# Save original diameters so each iteration is a fresh perturbation of the
# original tree, not a cumulative one.
orig_d_cm = copy(tree.segment_diameter_cm)

results = Tuple{Float64, Float64, Float64}[]  # (stenosis_pct, root_d_um, root_flow_mlmin)

for st in STENOSIS_PCTS
    tree.segment_diameter_cm .= orig_d_cm
    scale = 1.0 - st / 100.0
    for i in prox_idx
        tree.segment_diameter_cm[i] = orig_d_cm[i] * scale
    end
    t1 = time()
    hemo = compute_hemodynamics(tree;
        root_pressure=ROOT_P_PA,
        terminal_pressure=TERM_P_PA,
        hematocrit=HCT)
    t_hemo = time() - t1

    root_flow_m3s = 0.0
    for c in tree.children[tree.root_vertex]
        seg = tree.incoming_segment[c]
        seg != 0 && (root_flow_m3s += hemo.segment_flow[seg])
    end
    root_flow_mlmin = root_flow_m3s * 60e6
    stenosed_d_um = mean(orig_prox_d_um) * scale

    @printf("stenosis=%5.1f%%  proximal_d=%7.1f μm  root_flow=%7.2f mL/min  (hemo %.1fs)\n",
        st, stenosed_d_um, root_flow_mlmin, t_hemo)
    push!(results, (st, stenosed_d_um, root_flow_mlmin))
end

# Restore original diameters in-memory (polite, in case the tree is reused)
tree.segment_diameter_cm .= orig_d_cm

# Write CSV
mkpath(dirname(OUT_CSV))
open(OUT_CSV, "w") do io
    println(io, "stenosis_pct,proximal_diameter_um,root_flow_mlmin")
    for (st, d, q) in results
        @printf(io, "%.1f,%.4f,%.6f\n", st, d, q)
    end
end
println("-" ^ 80)
println("wrote $(OUT_CSV)")

# Also emit a quick ASCII chart for immediate feedback
println("\nASCII chart (root_flow vs stenosis%)")
max_q = maximum(r[3] for r in results)
for (st, _, q) in results
    bar = repeat("█", round(Int, 50 * q / max_q))
    @printf("  %4.0f%% | %s %.1f mL/min\n", st, bar, q)
end

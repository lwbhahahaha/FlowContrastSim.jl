#!/usr/bin/env julia
# LAD stenosis sweep: shrink the proximal LAD segments (label == "dias_lad1")
# by 0%, 10%, 20%, ..., 90% and report root flow at each step.
#
# The tree is loaded once; for each stenosis level we restore the original
# diameter vector, multiply the proximal indices, and re-solve Poiseuille +
# Pries + downstream cap_R.
#
# Stenosis is applied uniformly to every segment labeled "dias_lad1" — the
# proximal LAD chain (~12 mm). Diameter d_new = d_orig × (1 − stenosis/100).
#
# Usage:
#   julia --project=. scripts/lad_stenosis_sweep.jl \
#         <lad_csv_path> <config.toml> [output_csv]
#
# Example:
#   # baseline (at-rest LAD)
#   julia --project=. scripts/lad_stenosis_sweep.jl \
#         ../VascularTreeSim.jl/output_at_rest/lad_segments.csv \
#         configs/coronary_baseline.toml \
#         scripts/lad_stenosis_baseline.csv
#
#   # hyperemic (max-dilated LAD)
#   julia --project=. scripts/lad_stenosis_sweep.jl \
#         ../VascularTreeSim.jl/output/lad_segments.csv \
#         configs/coronary_hyperemic.toml \
#         scripts/lad_stenosis_hyperemic.csv

using FlowContrastSim
using Printf
using Statistics
import FlowContrastSim: load_tree, compute_hemodynamics, load_flow_config

length(ARGS) >= 2 || error("Usage: lad_stenosis_sweep.jl <lad_csv> <config.toml> [output_csv]")

const CSV_PATH    = ARGS[1]
const CONFIG_PATH = ARGS[2]
const cfg_base    = splitext(basename(CONFIG_PATH))[1]
const OUT_CSV     = length(ARGS) >= 3 ? ARGS[3] : "scripts/lad_stenosis_$(cfg_base).csv"

const PROXIMAL_LABEL = "dias_lad1"
const STENOSIS_PCTS = 0.0:10.0:90.0

cfg = load_flow_config(CONFIG_PATH)
mass_lad = get(cfg.territory_masses_g, "LAD", 0.0)
ROOT_P_PA = cfg.root_pressure_mmhg * 133.322
TERM_P_PA = cfg.terminal_pressure_mmhg * 133.322

println("LAD stenosis sweep")
println("  tree CSV: $(CSV_PATH)")
println("  config:   $(CONFIG_PATH)")
println("  output:   $(OUT_CSV)")
println("  proximal label: $(PROXIMAL_LABEL)")
println("  BCs: $(cfg.root_pressure_mmhg) → $(cfg.terminal_pressure_mmhg) mmHg, Hct=$(cfg.discharge_hematocrit)")
println("  cap_R=$(cfg.capillary_bed_R_per_100g_mmHgmin_ml) mmHg·min/mL/100g, LAD territory mass=$(mass_lad) g")
println("-" ^ 80)

print("loading LAD tree ... "); flush(stdout)
t0 = time()
tree = load_tree("LAD", CSV_PATH)
@printf("done in %.1fs, %d segs, %d verts\n", time()-t0, length(tree.segment_start), length(tree.vertices))

prox_idx = findall(==(PROXIMAL_LABEL), tree.segment_label)
if isempty(prox_idx)
    error("No segments with label '$(PROXIMAL_LABEL)' found in tree.")
end
orig_prox_d_um = [tree.segment_diameter_cm[i] * 1e4 for i in prox_idx]
@printf("  found %d proximal segments, diameter range [%.1f, %.1f] μm\n",
    length(prox_idx),
    minimum(orig_prox_d_um), maximum(orig_prox_d_um))

orig_d_cm = copy(tree.segment_diameter_cm)
results = Tuple{Float64, Float64, Float64, Float64}[]   # (stenosis%, prox_d_um, root_flow_mlmin, %flowing)

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
        hematocrit=cfg.discharge_hematocrit,
        capillary_bed_R_per_100g_mmHgmin_ml=cfg.capillary_bed_R_per_100g_mmHgmin_ml,
        territory_mass_g=mass_lad)
    t_hemo = time() - t1

    root_flow_m3s = 0.0
    for c in tree.children[tree.root_vertex]
        seg = tree.incoming_segment[c]
        seg != 0 && (root_flow_m3s += hemo.segment_flow[seg])
    end
    root_flow_mlmin = root_flow_m3s * 60e6
    stenosed_d_um = mean(orig_prox_d_um) * scale
    pct_flowing = 100.0 * count(>(1e-18), hemo.segment_flow) / length(hemo.segment_flow)

    @printf("stenosis=%5.1f%%  proximal_d=%7.1f μm  root_flow=%8.2f mL/min  flowing=%5.1f%%  (hemo %.1fs)\n",
        st, stenosed_d_um, root_flow_mlmin, pct_flowing, t_hemo)
    flush(stdout)
    push!(results, (st, stenosed_d_um, root_flow_mlmin, pct_flowing))
end

tree.segment_diameter_cm .= orig_d_cm
mkpath(dirname(OUT_CSV))
open(OUT_CSV, "w") do io
    println(io, "stenosis_pct,proximal_diameter_um,root_flow_mlmin,pct_flowing")
    for (st, d, q, fl) in results
        @printf(io, "%.1f,%.4f,%.6f,%.2f\n", st, d, q, fl)
    end
end
println("-" ^ 80)
println("wrote $(OUT_CSV)")

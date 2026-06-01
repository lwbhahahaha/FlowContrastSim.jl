#!/usr/bin/env julia
# End-to-end run on the 8 μm subdivided coronary tree using the Uniphase
# injection protocol + Taylor-Aris PDE path. Sparse mode at 50 μm to keep
# the per-segment matrix tractable; peak-only mode to skip the full time
# series (we only need the V2 snapshot for downstream basis_simulator).
#
# Outputs to FlowContrastSim.jl/output_8um_uniphase/.
#
# Usage:
#   julia --project=. scripts/run_8um_uniphase.jl

using FlowContrastSim
using Printf
using Dates
using LinearAlgebra
using Statistics

const TREE_DIR = joinpath(@__DIR__, "..", "..", "VascularTreeSim.jl", "output")
const OUT_DIR  = joinpath(@__DIR__, "..", "output_8um_uniphase")
const CONFIG   = joinpath(@__DIR__, "..", "configs", "coronary_hyperemic_uniphase.toml")
mkpath(OUT_DIR)

println("[8um] starting $(Dates.now())")
println("[8um] tree dir: $TREE_DIR")
println("[8um] config:   $CONFIG")
println("[8um] output:   $OUT_DIR")
flush(stdout)

cfg = load_flow_config(CONFIG)

# Build the AIF first so we can pick V2 acquisition time.
_, aif = protocol_to_aif(cfg.injection_protocol, cfg.patient_physiology;
                          dt=cfg.dt, t_max=cfg.t_end)
peak_idx = argmax(aif)
peak_t = (peak_idx - 1) * cfg.dt
# Clinical V2 is typically captured ~2 s after AIF peak (bolus tracking
# trigger + scan delay). We use peak + 2 s for the SVP V2 snapshot.
v2_t = peak_t + 2.0
@printf("[8um] AIF peak = %.3f mgI/mL at t = %.2f s, V2 snapshot at %.2f s\n",
        maximum(aif), peak_t, v2_t)
flush(stdout)

t_start = time()
result = run_flow_simulation(TREE_DIR, cfg;
                              output_dir=OUT_DIR,
                              peak_time_s=v2_t)
t_total = time() - t_start
@printf("[8um] total elapsed: %.1f min\n", t_total / 60.0)

# ── Dump per-tree summary ──
summary_path = joinpath(OUT_DIR, "summary.txt")
open(summary_path, "w") do io
    println(io, "FlowContrastSim 8 μm end-to-end run — $(Dates.now())")
    println(io, "=" ^ 70)
    println(io, "Tree dir: $TREE_DIR")
    println(io, "Config:   $CONFIG")
    println(io, "Protocol: $(typeof(result.protocol))")
    println(io, "  $(result.protocol)")
    println(io, "Physiology: $(result.physiology)")
    println(io, "AIF peak: $(round(maximum(aif), digits=3)) mgI/mL at t=$(round(peak_t,digits=2))s")
    println(io, "V2 snapshot at t=$(round(v2_t, digits=2))s")
    println(io, "Total elapsed: $(round(t_total / 60.0, digits=2)) min")
    println(io, "")
    for name in sort(collect(keys(result.trees)))
        tree = result.trees[name]
        cr   = result.contrast_results[name]
        hemo = result.hemo_results[name]
        n_seg = length(tree.segment_start)
        n_sim = isempty(cr.segment_ids) ? n_seg : length(cr.segment_ids)
        peak_v = maximum(cr.concentration)
        # Mean and max Taylor variance over reachable segs.
        σ²_finite = filter(isfinite, cr.arrival_variance_s2)
        σ_med = isempty(σ²_finite) ? NaN : sqrt(median(σ²_finite))
        σ_max = isempty(σ²_finite) ? NaN : sqrt(maximum(σ²_finite))
        arr_med = isempty(filter(isfinite, cr.arrival_s)) ? NaN :
                  median(filter(isfinite, cr.arrival_s))
        arr_max = isempty(filter(isfinite, cr.arrival_s)) ? NaN :
                  maximum(filter(isfinite, cr.arrival_s))
        root_flow_mL_min = 0.0
        for c in tree.children[tree.root_vertex]
            seg = tree.incoming_segment[c]
            seg != 0 && (root_flow_mL_min += hemo.segment_flow[seg])
        end
        root_flow_mL_min *= 60e6
        @printf(io, "[%s] %d segs total, %d simulated (≥%.0f μm)\n",
                name, n_seg, n_sim, cfg.contrast_min_diameter_um)
        @printf(io, "  root flow:           %.1f mL/min\n", root_flow_mL_min)
        @printf(io, "  arrival times (s):   median=%.2f, max=%.2f\n", arr_med, arr_max)
        @printf(io, "  Taylor σ (s):        median=%.3f, max=%.3f\n", σ_med, σ_max)
        @printf(io, "  V2 peak conc:        %.3f mgI/mL\n", peak_v)
    end
end
println("[8um] summary: $(summary_path)")
println("[8um] done $(Dates.now())")

#!/usr/bin/env julia
# Hemodynamics summary — skips contrast transport to avoid the O(N × T)
# concentration matrix (~45 GB per 30 M-segment tree at full resolution).
# Loads each tree once, then solves Poiseuille + Pries hemodynamics for any
# number of supplied config files (typically: baseline AND hyperemic), so
# both states can be verified from a single tree load.
#
# Usage:
#   julia --project=. scripts/natural_flow_summary.jl [tree_csv_dir] [config1.toml config2.toml ...]
#
# If no config files are passed, runs the legacy "no terminal R" mode and
# compares against hard-coded baseline/hyperemic targets.

using FlowContrastSim
using Printf
import FlowContrastSim: load_tree, compute_hemodynamics, load_flow_config

const TREE_DIR = length(ARGS) >= 1 ? ARGS[1] :
    "../VascularTreeSim.jl/output"
const CONFIG_PATHS = length(ARGS) >= 2 ? ARGS[2:end] : String[]

# Legacy targets (used only when no configs supplied)
const LEGACY_BASELINE  = Dict("LAD" => 50.0, "LCX" => 30.0, "RCA" => 50.0)
const LEGACY_HYPEREMIC = Dict("LAD" => 200.0, "LCX" => 150.0, "RCA" => 190.0)

# Load configs (or build a single "no terminal R" legacy config)
configs = NamedTuple[]
if isempty(CONFIG_PATHS)
    push!(configs, (
        name = "no_terminal_R",
        root_p_pa = 100.0 * 133.322,
        term_p_pa = 15.0 * 133.322,
        hct = 0.45,
        cap_R = 0.0,
        masses = Dict{String,Float64}(),
        targets = Dict{String,Float64}(),
    ))
else
    for path in CONFIG_PATHS
        cfg = load_flow_config(path)
        push!(configs, (
            name = basename(path),
            root_p_pa = cfg.root_pressure_mmhg * 133.322,
            term_p_pa = cfg.terminal_pressure_mmhg * 133.322,
            hct = cfg.discharge_hematocrit,
            cap_R = cfg.capillary_bed_R_per_100g_mmHgmin_ml,
            masses = cfg.territory_masses_g,
            targets = cfg.target_flows_ml_min,
        ))
    end
end

println("Hemodynamics summary  tree_dir=$(TREE_DIR)")
for c in configs
    println("  config=$(c.name)  root=$(c.root_p_pa/133.322) term=$(c.term_p_pa/133.322) Hct=$(c.hct)  cap_R=$(c.cap_R) mmHg·min/mL/100g")
end
println("-" ^ 80)

# results[config_name] => Vector of NamedTuples
results = Dict{String, Vector{NamedTuple}}()
for c in configs
    results[c.name] = NamedTuple[]
end

for name in ("LAD", "LCX", "RCA")
    csv = joinpath(TREE_DIR, "$(lowercase(name))_segments.csv")
    if !isfile(csv)
        println("[$(name)] CSV missing: $(csv)"); continue
    end
    t0 = time()
    print("[$(name)] loading ..."); flush(stdout)
    tree = load_tree(name, csv)
    t_load = time() - t0
    println(" loaded ($(round(t_load; digits=1))s, $(length(tree.segment_start)) segs)")
    flush(stdout)

    for c in configs
        t1 = time()
        mass = get(c.masses, name, 0.0)
        hemo = compute_hemodynamics(tree;
            root_pressure=c.root_p_pa,
            terminal_pressure=c.term_p_pa,
            hematocrit=c.hct,
            capillary_bed_R_per_100g_mmHgmin_ml=c.cap_R,
            territory_mass_g=mass)
        t_hemo = time() - t1

        root_flow_m3s = 0.0
        for ch in tree.children[tree.root_vertex]
            seg = tree.incoming_segment[ch]
            seg != 0 && (root_flow_m3s += hemo.segment_flow[seg])
        end
        root_flow_mlmin = root_flow_m3s * 60e6

        n_flow = count(>(1e-18), hemo.segment_flow)
        pct_flow = round(100 * n_flow / length(hemo.segment_flow); digits=1)

        target = get(c.targets, name, 0.0)
        ratio_str = target > 0 ? @sprintf("Q/target=%.2f", root_flow_mlmin/target) : "no target"
        println("  [$(c.name)] hemo $(round(t_hemo; digits=1))s | flow=$(round(root_flow_mlmin; digits=1)) mL/min | $(pct_flow)% segs flowing | $(ratio_str)")
        flush(stdout)

        push!(results[c.name], (name=name, segs=length(tree.segment_start),
                                flow=root_flow_mlmin, target=target,
                                pct_flow=pct_flow, t_hemo=t_hemo))
        hemo = nothing
        GC.gc()
    end
    # Release tree before next iter
    tree = nothing
    GC.gc()
end

println("-" ^ 80)
for c in configs
    println("SUMMARY ($(c.name))")
    @printf("%-5s %10s %12s %10s %10s\n", "tree", "segs", "flow_mlmin", "target", "%flowing")
    total = 0.0
    for r in results[c.name]
        @printf("%-5s %10d %12.1f %10.1f %10.1f\n", r.name, r.segs, r.flow, r.target, r.pct_flow)
        total += r.flow
    end
    @printf("%-5s %10s %12.1f\n", "SUM", "", total)
    println()
end

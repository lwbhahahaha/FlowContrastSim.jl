#!/usr/bin/env julia
# Minimal natural-physics hemodynamics summary — skips contrast transport to
# avoid the O(N × T) concentration matrix (~45 GB per 30 M-segment tree at
# full resolution). Just loads each tree, solves Poiseuille + Pries, reports
# root flow vs baseline/hyperemic targets.
#
# Usage:
#   julia --project=. scripts/natural_flow_summary.jl [tree_csv_dir]

using FlowContrastSim
using Printf
import FlowContrastSim: load_tree, compute_hemodynamics

const TREE_DIR = length(ARGS) >= 1 ? ARGS[1] :
    "../VascularTreeSim.jl/output"

# Baseline (at rest) targets. Hyperemic = baseline × 1.6^3 ≈ × 4.096.
const BASELINE = Dict("LAD" => 59.08, "LCX" => 28.32, "RCA" => 52.25)
const HYPEREMIC = Dict("LAD" => 242.0, "LCX" => 116.0, "RCA" => 214.0)

const ROOT_P_PA = 100.0 * 133.322
const TERM_P_PA = 15.0 * 133.322
const HCT = 0.45

println("Natural-physics hemo summary  tree_dir=$(TREE_DIR)")
println("BCs: root=$(ROOT_P_PA/133.322) mmHg, term=$(TERM_P_PA/133.322) mmHg, Hct=$(HCT)")
println("-" ^ 80)

results = NamedTuple[]
for name in ("LAD", "LCX", "RCA")
    csv = joinpath(TREE_DIR, "$(lowercase(name))_segments.csv")
    if !isfile(csv)
        println("[$(name)] CSV missing: $(csv)"); continue
    end
    t0 = time()
    print("[$(name)] loading ..."); flush(stdout)
    tree = load_tree(name, csv)
    t_load = time() - t0
    print(" loaded ($(round(t_load; digits=1))s, $(length(tree.segment_start)) segs) ... ")
    flush(stdout)

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

    n_flow = count(>(1e-18), hemo.segment_flow)
    pct_flow = round(100 * n_flow / length(hemo.segment_flow); digits=1)

    base = BASELINE[name]; hyper = HYPEREMIC[name]
    println("hemo $(round(t_hemo; digits=1))s | flow=$(round(root_flow_mlmin; digits=1)) mL/min | $(pct_flow)% segs flowing | base×$(round(root_flow_mlmin/base; digits=2)) | hyper×$(round(root_flow_mlmin/hyper; digits=2))")
    push!(results, (name=name, segs=length(tree.segment_start),
                    flow=root_flow_mlmin, base=base, hyper=hyper,
                    pct_flow=pct_flow, t_load=t_load, t_hemo=t_hemo))
    # Release tree + hemo before next iter
    tree = nothing; hemo = nothing
    GC.gc()
end

println("-" ^ 80)
println("SUMMARY (natural, no terminal-resistance calibration)")
@printf("%-5s %10s %12s %10s %10s %10s\n", "tree", "segs", "flow_mlmin", "baseline", "hyperemic", "%flowing")
for r in results
    @printf("%-5s %10d %12.1f %10.1f %10.1f %10.1f\n", r.name, r.segs, r.flow, r.base, r.hyper, r.pct_flow)
end

#!/usr/bin/env julia
# Diagnostic: characterize the LCX tree's dead-flow segments.
#
# Reports for each tree:
#   - Total segs, flowing fraction (above 1e-18 m³/s) and "near-dead" fractions
#   - For DEAD segments (flow < 1e-22 m³/s): histograms by
#       * diameter (μm log10 bins)
#       * length (mm log10 bins)
#       * label (grown / subdivided / xcat)
#       * BFS depth from root
#   - Top-10 "bottleneck" segments: segments whose flow / parent_flow ratio is
#     extreme (i.e., a parallel split heavily favors the sibling). These are
#     the gateways to dead subtrees.
#
# Usage:
#   julia --project=. scripts/lcx_dead_flow_diag.jl [tree_csv_dir] [tree_name]

using FlowContrastSim
using Printf
using Statistics
import FlowContrastSim: load_tree, compute_hemodynamics, _topo_order

const TREE_DIR  = length(ARGS) >= 1 ? ARGS[1] : "../VascularTreeSim.jl/output"
const TREE_NAME = length(ARGS) >= 2 ? uppercase(ARGS[2]) : "LCX"

const ROOT_P_PA = 100.0 * 133.322
const TERM_P_PA = 15.0 * 133.322
const HCT = 0.45

const FLOW_DEAD_THRESH = 1e-22  # essentially zero (m³/s)
const FLOW_FLOW_THRESH = 1e-18  # threshold used by compute_hemodynamics "flowing" stat

function bin_log10(x::Float64, edges::Vector{Float64})
    x <= 0 && return 1
    lx = log10(x)
    idx = searchsortedfirst(edges, lx)
    return clamp(idx, 1, length(edges) + 1)
end

csv = joinpath(TREE_DIR, "$(lowercase(TREE_NAME))_segments.csv")
isfile(csv) || error("Tree CSV missing: $csv")

@printf("Loading %s tree from %s ...\n", TREE_NAME, csv); flush(stdout)
t0 = time()
tree = load_tree(TREE_NAME, csv)
@printf("  loaded in %.1fs: %d segs, %d verts\n", time()-t0,
        length(tree.segment_start), length(tree.vertices)); flush(stdout)

@printf("Computing hemodynamics (no cap R) ...\n"); flush(stdout)
t1 = time()
hemo = compute_hemodynamics(tree;
    root_pressure=ROOT_P_PA, terminal_pressure=TERM_P_PA, hematocrit=HCT)
@printf("  hemo done in %.1fs\n", time()-t1); flush(stdout)

nseg = length(tree.segment_start)
flow = hemo.segment_flow

# ── Summary stats ──
n_dead   = count(<(FLOW_DEAD_THRESH), flow)
n_flow18 = count(<(FLOW_FLOW_THRESH), flow)
n_real   = nseg - n_flow18
println()
println("-" ^ 80)
@printf("Flow summary for %s tree:\n", TREE_NAME)
@printf("  total segs:                          %d\n", nseg)
@printf("  flowing (flow >= 1e-18 m³/s):        %d  (%.2f%%)\n", n_real, 100*n_real/nseg)
@printf("  below 1e-18 m³/s (\"dead\" by hemo):   %d  (%.2f%%)\n", n_flow18, 100*n_flow18/nseg)
@printf("  below 1e-22 m³/s (essentially zero): %d  (%.2f%%)\n", n_dead, 100*n_dead/nseg)
println()

# ── Dead-segment diagnostics ──
order = _topo_order(tree)
depth = Dict{Int, Int}()
# Robust depth: tree may have merges where incoming_segment[sv] doesn't match
# the BFS parent. Use get() with a sentinel; segments whose parent isn't yet
# recorded (or which are themselves orphans) are flagged depth = -1.
for s in order
    sv = tree.segment_start[s]
    parent_seg = tree.incoming_segment[sv]
    if parent_seg == 0
        depth[s] = 0
    else
        pd = get(depth, parent_seg, -1)
        depth[s] = pd < 0 ? -1 : pd + 1
    end
end

dead_idx = findall(<(FLOW_FLOW_THRESH), flow)
isempty(dead_idx) && (println("No dead segments found — done."); exit())

# Distribution of dead segments
diameters_um = [tree.segment_diameter_cm[i] * 1e4 for i in dead_idx]
lengths_mm = Float64[]
for i in dead_idx
    a = tree.vertices[tree.segment_start[i]]
    b = tree.vertices[tree.segment_end[i]]
    push!(lengths_mm, 10.0 * sqrt(sum((a .- b).^2)))
end
depths_dead = [get(depth, i, -1) for i in dead_idx]
labels_dead = [tree.segment_label[i] for i in dead_idx]

println("-" ^ 80)
println("Dead-segment diameter histogram (log10 μm bins):")
edges_d = collect(0.5:0.25:4.0)
hist_d = zeros(Int, length(edges_d) + 1)
for d in diameters_um
    hist_d[bin_log10(d, edges_d)] += 1
end
@printf("  %-15s %12s %8s\n", "log10_um", "count", "pct")
for i in 1:length(edges_d)
    a = i == 1 ? -Inf : edges_d[i-1]
    b = edges_d[i]
    @printf("  [%5.2f,%5.2f)   %12d %7.1f%%\n", a, b, hist_d[i], 100*hist_d[i]/n_flow18)
end
@printf("  [>=%4.2f]         %12d %7.1f%%\n", edges_d[end], hist_d[end], 100*hist_d[end]/n_flow18)

println()
println("Dead-segment length histogram (log10 mm bins):")
edges_l = collect(-4.0:0.5:1.0)
hist_l = zeros(Int, length(edges_l) + 1)
for L in lengths_mm
    hist_l[bin_log10(L, edges_l)] += 1
end
@printf("  %-15s %12s %8s\n", "log10_mm", "count", "pct")
for i in 1:length(edges_l)
    a = i == 1 ? -Inf : edges_l[i-1]
    b = edges_l[i]
    @printf("  [%5.2f,%5.2f)   %12d %7.1f%%\n", a, b, hist_l[i], 100*hist_l[i]/n_flow18)
end
@printf("  [>=%4.2f]         %12d %7.1f%%\n", edges_l[end], hist_l[end], 100*hist_l[end]/n_flow18)

println()
println("Dead-segment label distribution:")
label_counts = Dict{String, Int}()
for lbl in labels_dead
    label_counts[lbl] = get(label_counts, lbl, 0) + 1
end
for (lbl, c) in sort(collect(label_counts), by=x -> -x[2])[1:min(20, length(label_counts))]
    @printf("  %-30s %12d  %5.1f%%\n", lbl, c, 100*c/n_flow18)
end

println()
println("Dead-segment depth (BFS from root) distribution:")
edges_dp = collect(0:50:1500)
hist_dp = zeros(Int, length(edges_dp) + 1)
for d in depths_dead
    d < 0 && continue
    idx = searchsortedfirst(edges_dp, d)
    idx = clamp(idx, 1, length(edges_dp) + 1)
    hist_dp[idx] += 1
end
for i in 1:length(edges_dp)
    hist_dp[i] == 0 && continue
    a = i == 1 ? 0 : edges_dp[i-1]
    b = edges_dp[i]
    @printf("  [%4d,%4d)   %12d %7.1f%%\n", a, b, hist_dp[i], 100*hist_dp[i]/n_flow18)
end
if hist_dp[end] > 0
    @printf("  [>=%4d]       %12d %7.1f%%\n", edges_dp[end], hist_dp[end], 100*hist_dp[end]/n_flow18)
end

# ── Bottleneck identification ──
# A "bottleneck" is a segment where flow drops dramatically relative to parent.
# We look at child/parent flow ratios; small ratio = this branch loses out.
println()
println("-" ^ 80)
println("Top-30 \"bottleneck\" segments (smallest flow / parent_flow ratio, parent flow > 1e-12):")
bottleneck = Tuple{Float64, Int, Int}[]   # (ratio, seg, parent_seg)
for s in 1:nseg
    sv = tree.segment_start[s]
    parent_seg = tree.incoming_segment[sv]
    parent_seg == 0 && continue
    pf = flow[parent_seg]
    pf <= 1e-12 && continue
    cf = flow[s]
    ratio = cf / pf
    push!(bottleneck, (ratio, s, parent_seg))
end
sort!(bottleneck, by=first)
n_show = min(30, length(bottleneck))
@printf("  %-8s %-12s %-12s %-12s %-12s %-12s %-20s\n",
        "ratio", "seg_id", "depth", "d_um", "L_mm", "p_flow_mlmin", "label")
for k in 1:n_show
    r, s, ps = bottleneck[k]
    d_um = tree.segment_diameter_cm[s] * 1e4
    a = tree.vertices[tree.segment_start[s]]
    b = tree.vertices[tree.segment_end[s]]
    L_mm = 10.0 * sqrt(sum((a .- b).^2))
    pf_ml = flow[ps] * 60e6
    @printf("  %.2e  %-12d %-12d %-12.2f %-12.4f %-12.2e %-20s\n",
            r, s, get(depth, s, -1), d_um, L_mm, pf_ml, tree.segment_label[s])
end

# How many child/parent ratios are very small? indicates many "lossy" bifurcations
n_skew_1k = count(b -> b[1] < 1e-3, bottleneck)
n_skew_1M = count(b -> b[1] < 1e-6, bottleneck)
n_skew_1B = count(b -> b[1] < 1e-9, bottleneck)
println()
@printf("Bifurcation skew (ratio = child_flow / parent_flow):\n")
@printf("  ratio < 1e-3 (>1000x split): %d  (%.2f%%)\n", n_skew_1k, 100*n_skew_1k/length(bottleneck))
@printf("  ratio < 1e-6:                %d  (%.2f%%)\n", n_skew_1M, 100*n_skew_1M/length(bottleneck))
@printf("  ratio < 1e-9:                %d  (%.2f%%)\n", n_skew_1B, 100*n_skew_1B/length(bottleneck))

println()
println("done.")

#!/usr/bin/env julia
# Per-tree conductance diagnostic. Answers: why does LCX/RCA natural flow
# exceed target by 2-2.5× while LAD is only 1.2× over? Compares:
#   - root segment geometry (d, L, single-seg R)
#   - path statistics (avg/min depth, avg/min root→terminal R sum)
#   - diameter histogram across 8 bins from smallest to root
#   - n_terminals, n_segments, avg fanout
#
# No calibration, no contrast — just topological/geometric analysis.
# Reads ../VascularTreeSim.jl/output/<tree>_segments.csv directly.

using FlowContrastSim
using Printf
using Statistics
import FlowContrastSim: load_tree, _compute_resistance, _topo_order

const TREE_DIR = length(ARGS) >= 1 ? ARGS[1] :
    "../VascularTreeSim.jl/output"
const HCT = 0.45

println("Tree diagnostic  dir=$(TREE_DIR)  Hct=$(HCT)")
println("-" ^ 100)

for name in ("LAD", "LCX", "RCA")
    csv = joinpath(TREE_DIR, "$(lowercase(name))_segments.csv")
    isfile(csv) || (println("[$(name)] MISSING"); continue)

    print("[$(name)] loading... "); flush(stdout)
    tree = load_tree(name, csv)
    nseg = length(tree.segment_start)
    nvert = length(tree.vertices)
    println("$(nseg) segs, $(nvert) verts")

    # Resistance per segment (Pries)
    R = _compute_resistance(tree; hematocrit=HCT)

    # Root segment geometry
    root_v = tree.root_vertex
    root_children = tree.children[root_v]
    root_segs = Int[]
    for c in root_children
        s = tree.incoming_segment[c]
        s != 0 && push!(root_segs, s)
    end
    println("  root vertex has $(length(root_children)) children, $(length(root_segs)) root segments")
    for rs in root_segs
        sv = tree.segment_start[rs]; ev = tree.segment_end[rs]
        a = tree.vertices[sv]; b = tree.vertices[ev]
        seg_L_m = sqrt(sum((a .- b).^2)) * 0.01   # cm → m
        seg_d_m = tree.segment_diameter_cm[rs] * 0.01
        println("    root_seg#$(rs): d=$(round(seg_d_m*1e6; digits=1))μm L=$(round(seg_L_m*1e3; digits=2))mm R=$(round(R[rs]; sigdigits=3)) Pa·s/m³")
    end

    # Diameter histogram in μm
    d_um = tree.segment_diameter_cm .* 1e4
    d_min = minimum(d_um); d_max = maximum(d_um)
    @printf("  diameter range: [%.1f, %.1f] μm (mean=%.1f, median=%.1f)\n",
        d_min, d_max, mean(d_um), median(d_um))

    # Histogram: log-bins from min to max
    edges = exp.(range(log(d_min), log(d_max); length=9))
    counts = zeros(Int, 8)
    for d in d_um
        for i in 1:8
            if d <= edges[i+1]
                counts[i] += 1; break
            end
        end
    end
    println("  diameter log-bins:")
    for i in 1:8
        bar_len = round(Int, 60 * counts[i] / nseg)
        @printf("    %6.1f - %6.1f μm : %10d (%5.1f%%) %s\n",
            edges[i], edges[i+1], counts[i], 100*counts[i]/nseg, "█"^bar_len)
    end

    # Path to root: for each terminal, count hops and sum R along path
    parent_seg = zeros(Int, nseg)
    for s in 1:nseg
        sv = tree.segment_start[s]
        parent_seg[s] = tree.incoming_segment[sv]
    end

    # Terminals = segments whose end vertex has no outgoing children
    terminal_segs = Int[]
    for s in 1:nseg
        ev = tree.segment_end[s]
        if isempty(tree.children[ev])
            push!(terminal_segs, s)
        end
    end
    println("  terminals: $(length(terminal_segs))")

    # For efficiency, sample up to 10_000 terminals
    sample_terms = length(terminal_segs) > 10_000 ?
        terminal_segs[rand(1:length(terminal_segs), 10_000)] : terminal_segs
    path_hops = Int[]
    path_R = Float64[]
    for t in sample_terms
        hops = 0
        Rsum = 0.0
        s = t
        while s != 0 && hops < 1_000_000
            hops += 1
            Rsum += R[s]
            s = parent_seg[s]
        end
        push!(path_hops, hops)
        push!(path_R, Rsum)
    end
    @printf("  path hops  (terminals→root, sampled %d): mean=%.1f median=%.0f min=%d max=%d\n",
        length(sample_terms), mean(path_hops), median(path_hops),
        minimum(path_hops), maximum(path_hops))
    # R in Pa·s/m³; express total path R more intuitively as ΔP/Q at 85 mmHg
    # If a single terminal's path R determines its flow (crude, no merges),
    # Q_term = 85mmHg / Rsum. Report as mL/min-equivalent.
    dP_pa = 85.0 * 133.322
    q_term_mlmin = (dP_pa ./ path_R) .* 60e6
    @printf("  path R-sum (per terminal, if isolated): mean=%.2e min=%.2e max=%.2e Pa·s/m³\n",
        mean(path_R), minimum(path_R), maximum(path_R))
    @printf("  per-terminal Q if isolated @ 85mmHg: mean=%.3f min=%.3f max=%.3f mL/min\n",
        mean(q_term_mlmin), minimum(q_term_mlmin), maximum(q_term_mlmin))

    tree = nothing; R = nothing; GC.gc()
    println("-" ^ 100)
end

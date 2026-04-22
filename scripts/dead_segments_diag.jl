#!/usr/bin/env julia
# Dead-segment forensic: compare two connectivity measures per tree:
#   (a) parent_segment_id chain (raw CSV) — how many segs reach a root by
#       walking parent_id
#   (b) FlowTree BFS from root_vertex — how many segs are in _topo_order
#
# If (a) > (b), load_tree's vertex-coincidence logic is losing segments
# (because multiple segs share an end-vertex, so `incoming_segment[ev]`
# gets overwritten).
#
# Also reports how many distinct segments end at each vertex (merge
# factor). A healthy tree has merge factor ≤ 1 per vertex.

using FlowContrastSim
using Printf
import FlowContrastSim: load_tree, _topo_order

const TREE_DIR = length(ARGS) >= 1 ? ARGS[1] :
    "../VascularTreeSim.jl/output"

function reachable_via_parentid(csv_path::String)
    # Minimal CSV parse: segment_id, parent_segment_id columns only.
    io = open(csv_path, "r")
    try
        header = split(readline(io), ',')
        sid_col = findfirst(==("segment_id"), header)
        pid_col = findfirst(==("parent_segment_id"), header)
        sid_col === nothing && return (0, 0, 0, Int[])
        pid_col === nothing && return (0, 0, 0, Int[])

        parent_of = Dict{Int, Int}()
        all_ids = Int[]
        roots = Int[]
        while !eof(io)
            line = readline(io); isempty(line) && continue
            fields = split(line, ',')
            sid = parse(Int, fields[sid_col])
            pid = isempty(fields[pid_col]) ? 0 : parse(Int, fields[pid_col])
            parent_of[sid] = pid
            push!(all_ids, sid)
            pid <= 0 && push!(roots, sid)
        end

        # BFS from roots via reverse parent_of
        children_of = Dict{Int, Vector{Int}}()
        for (sid, pid) in parent_of
            if pid > 0
                push!(get!(children_of, pid, Int[]), sid)
            end
        end
        reached = Set{Int}()
        queue = copy(roots)
        head = 1
        while head <= length(queue)
            s = queue[head]; head += 1
            s in reached && continue
            push!(reached, s)
            for c in get(children_of, s, Int[])
                c in reached || push!(queue, c)
            end
        end
        return (length(all_ids), length(reached), length(roots), roots)
    finally
        close(io)
    end
end

for name in ("LAD", "LCX", "RCA")
    csv = joinpath(TREE_DIR, "$(lowercase(name))_segments.csv")
    isfile(csv) || (println("[$(name)] MISSING"); continue)

    println("=" ^ 80)
    println("[$(name)] CSV: $(csv)")

    # (a) parent_segment_id reachability
    print("  (a) walking parent_segment_id chain... "); flush(stdout)
    t0 = time()
    (nseg_csv, n_reached_csv, nroots_csv, roots_csv) = reachable_via_parentid(csv)
    @printf("done in %.1fs\n", time()-t0)
    @printf("      CSV segs: %d, reached via parent_id: %d (%.2f%%), roots: %d\n",
        nseg_csv, n_reached_csv, 100*n_reached_csv/nseg_csv, nroots_csv)

    # (b) FlowTree BFS
    print("  (b) loading FlowTree + BFS... "); flush(stdout)
    t1 = time()
    tree = load_tree(name, csv)
    order = _topo_order(tree)
    @printf("done in %.1fs\n", time()-t1)
    nseg_tree = length(tree.segment_start)
    @printf("      FlowTree segs: %d, reached via BFS: %d (%.2f%%)\n",
        nseg_tree, length(order), 100*length(order)/nseg_tree)

    # Merge analysis — how many distinct segs end at the same vertex?
    merge_count = Dict{Int, Int}()
    for i in 1:nseg_tree
        ev = tree.segment_end[i]
        merge_count[ev] = get(merge_count, ev, 0) + 1
    end
    merge_histo = Dict{Int, Int}()
    for (_, c) in merge_count
        merge_histo[c] = get(merge_histo, c, 0) + 1
    end
    println("  vertex merge histogram (how many segs end at each vertex):")
    for k in sort(collect(keys(merge_histo)))
        @printf("      %3d segs/vertex : %d vertices\n", k, merge_histo[k])
    end

    # Start-vertex sharing: how many distinct segs START at the same vertex?
    start_count = Dict{Int, Int}()
    for i in 1:nseg_tree
        sv = tree.segment_start[i]
        start_count[sv] = get(start_count, sv, 0) + 1
    end
    start_histo = Dict{Int, Int}()
    for (_, c) in start_count
        start_histo[c] = get(start_histo, c, 0) + 1
    end
    println("  vertex bifurcation histogram (how many segs START at each vertex):")
    for k in sort(collect(keys(start_histo)))
        if k >= 3 || get(start_histo, k, 0) > 1_000_000
            @printf("      %3d segs/vertex : %d vertices\n", k, start_histo[k])
        end
    end
    @printf("      (1-2 segs/vertex is normal tree bifurcation; >=3 suggests merges)\n")

    tree = nothing; GC.gc()
    println()
end

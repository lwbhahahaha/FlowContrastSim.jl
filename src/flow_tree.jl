# ── FlowTree: read-only vascular tree loaded from CSV ──
# This is the interface between tree-growing and flow simulation.
# CSV is the contract: tree-growing produces it, liquid-flowing consumes it.

struct FlowTree
    name::String
    vertices::Vector{SVector{3,Float64}}
    parent_vertex::Vector{Int}
    children::Vector{Vector{Int}}
    incoming_segment::Vector{Int}
    segment_start::Vector{Int}
    segment_end::Vector{Int}
    segment_diameter_cm::Vector{Float64}
    segment_label::Vector{String}
    root_vertex::Int
end

"""
    load_tree(name, csv_path) -> FlowTree

Load a vascular tree from a CSV file produced by the tree-growing repo.

CSV columns: branch, segment_id, parent_segment_id, x1_cm, y1_cm, z1_cm,
             x2_cm, y2_cm, z2_cm, xmid_cm, ymid_cm, zmid_cm,
             length_mm, diameter_um, label

Topology source of truth: `parent_segment_id`. The grower writes a
deterministic chain (each segment's parent_segment_id points to the segment
whose end-vertex was this segment's start-vertex), so we rebuild the tree
purely from that chain rather than coordinate coincidence.

Why this matters: when segments are densely packed (~110 M for an LCX subdivided
to 6 μm), unrelated end-points can round to the same 1e-8 cm grid (10 nm), and
a vertex-coincidence dedup collapses them into one vertex. With only one
`incoming_segment[v]` slot per vertex, every collision orphans one of the
colliding segments + its entire downstream subtree. This created the
LCX "21% dead flow" we observed: a few thousand grown segs were orphaned and
took ~23 M subdivided descendants with them.

Vertices are still stored to provide segment endpoint coordinates for the
hemodynamics length calculation; only the topology is parent-id-driven.
The CSV's `parent_segment_id == 0` (or empty / negative) marks roots.
"""
function load_tree(name::String, csv_path::String)
    io = open(csv_path, "r")
    try
        header_line = readline(io)
        hdr = String.(split(header_line, ','))
        col(n) = findfirst(==(n), hdr)
        x1i = col("x1_cm"); y1i = col("y1_cm"); z1i = col("z1_cm")
        x2i = col("x2_cm"); y2i = col("y2_cm"); z2i = col("z2_cm")
        di  = col("diameter_um")
        li  = col("label")
        sidi = col("segment_id")
        psidi = col("parent_segment_id")
        any(x -> x === nothing, (x1i, y1i, z1i, x2i, y2i, z2i, di, li)) &&
            error("CSV $(csv_path) is missing required columns; got header: $(hdr)")
        has_parent_seg = sidi !== nothing && psidi !== nothing

        # Per-segment endpoint coords + scalar fields.
        sp_arr = SVector{3,Float64}[]
        ep_arr = SVector{3,Float64}[]
        seg_diam_cm = Float64[]
        seg_label   = String[]
        seg_id_csv  = Int[]
        parent_seg_csv = Int[]
        sizehint!(sp_arr, 4_000_000)
        sizehint!(ep_arr, 4_000_000)
        sizehint!(seg_diam_cm, 4_000_000)
        sizehint!(seg_label, 4_000_000)
        has_parent_seg && (sizehint!(seg_id_csv, 4_000_000); sizehint!(parent_seg_csv, 4_000_000))

        fields = Vector{SubString{String}}(undef, length(hdr))

        @inbounds while !eof(io)
            line = readline(io)
            isempty(line) && continue
            nf = 0
            last = 1
            lline = length(line)
            for idx in 1:lline
                @inbounds c = line[idx]
                if c == ','
                    nf += 1
                    fields[nf] = SubString(line, last, idx-1)
                    last = idx + 1
                end
            end
            nf += 1
            fields[nf] = SubString(line, last, lline)

            push!(sp_arr, SVector(parse(Float64, fields[x1i]),
                                  parse(Float64, fields[y1i]),
                                  parse(Float64, fields[z1i])))
            push!(ep_arr, SVector(parse(Float64, fields[x2i]),
                                  parse(Float64, fields[y2i]),
                                  parse(Float64, fields[z2i])))
            push!(seg_diam_cm, parse(Float64, fields[di]) * 1e-4)
            push!(seg_label,   String(fields[li]))
            if has_parent_seg
                push!(seg_id_csv, parse(Int, fields[sidi]))
                pstr = fields[psidi]
                pv = isempty(pstr) ? -1 : parse(Int, pstr)
                push!(parent_seg_csv, pv)
            end
        end

        nseg = length(sp_arr)
        nseg == 0 && error("CSV $(csv_path) contains no segment rows")

        if has_parent_seg
            return _build_flowtree_from_parent_chain(name, sp_arr, ep_arr,
                seg_diam_cm, seg_label, seg_id_csv, parent_seg_csv, csv_path)
        else
            return _build_flowtree_legacy(name, sp_arr, ep_arr,
                seg_diam_cm, seg_label, csv_path)
        end
    finally
        close(io)
    end
end

# Parent-id-driven topology reconstruction. Each segment is given its own
# unique end vertex; its start vertex is its parent's end vertex (or a fresh
# root vertex). This eliminates vertex-coincidence overwrites of
# `incoming_segment[v]`.
function _build_flowtree_from_parent_chain(name::String,
        sp_arr::Vector{SVector{3,Float64}}, ep_arr::Vector{SVector{3,Float64}},
        seg_diam_cm::Vector{Float64}, seg_label::Vector{String},
        seg_id_csv::Vector{Int}, parent_seg_csv::Vector{Int},
        csv_path::String)
    nseg = length(sp_arr)

    seg_id_to_idx = Dict{Int, Int}()
    sizehint!(seg_id_to_idx, nseg)
    for i in 1:nseg
        seg_id_to_idx[seg_id_csv[i]] = i
    end

    # Vertex allocation:
    #   - 1 vertex per segment for its end (unique).
    #   - For roots / orphans (parent not in the CSV), allocate a fresh start
    #     vertex too. Non-root segments share their parent's end vertex.
    vertices = SVector{3,Float64}[]
    sizehint!(vertices, nseg + 16)
    seg_end_vid = zeros(Int, nseg)
    seg_start_vid = zeros(Int, nseg)

    # Pass 1: end vertices.
    for i in 1:nseg
        push!(vertices, ep_arr[i])
        seg_end_vid[i] = length(vertices)
    end

    # Pass 2: start vertices via parent chain.
    n_orphans = 0
    root_segs = Int[]
    for i in 1:nseg
        pid = parent_seg_csv[i]
        if pid <= 0
            push!(vertices, sp_arr[i])
            seg_start_vid[i] = length(vertices)
            push!(root_segs, i)
        else
            parent_idx = get(seg_id_to_idx, pid, 0)
            if parent_idx > 0
                seg_start_vid[i] = seg_end_vid[parent_idx]
            else
                # Parent id referenced but not in CSV — treat as detached root.
                push!(vertices, sp_arr[i])
                seg_start_vid[i] = length(vertices)
                n_orphans += 1
            end
        end
    end

    nv = length(vertices)
    parent_vertex    = zeros(Int, nv)
    children         = [Int[] for _ in 1:nv]
    incoming_segment = zeros(Int, nv)

    for i in 1:nseg
        sv = seg_start_vid[i]
        ev = seg_end_vid[i]
        parent_vertex[ev] = sv
        incoming_segment[ev] = i  # unique by construction
        push!(children[sv], ev)
    end

    if isempty(root_segs)
        error("Could not locate a root segment in $(csv_path): no parent_segment_id <= 0")
    end
    # Prefer the first root whose start vertex is genuinely inflow-free.
    root_vertex = 0
    for r in root_segs
        rv = seg_start_vid[r]
        if incoming_segment[rv] == 0
            root_vertex = rv
            break
        end
    end
    root_vertex == 0 && (root_vertex = seg_start_vid[root_segs[1]])

    if n_orphans > 0
        @warn "load_tree: $(n_orphans) segments reference a parent_segment_id not present in CSV; treating them as detached roots." csv=csv_path
    end

    return FlowTree(name, vertices, parent_vertex, children, incoming_segment,
                    seg_start_vid, seg_end_vid, seg_diam_cm, seg_label, root_vertex)
end

# Fallback for legacy CSVs without segment_id / parent_segment_id. Same
# vertex-coincidence dedup as the original implementation; only used when the
# CSV doesn't supply the parent chain.
function _build_flowtree_legacy(name::String,
        sp_arr::Vector{SVector{3,Float64}}, ep_arr::Vector{SVector{3,Float64}},
        seg_diam_cm::Vector{Float64}, seg_label::Vector{String},
        csv_path::String)
    nseg = length(sp_arr)

    vertex_map = Dict{SVector{3,Float64}, Int}()
    vertices = SVector{3,Float64}[]
    sizehint!(vertex_map, 4_000_000)
    sizehint!(vertices, 4_000_000)
    function get_or_add!(pt::SVector{3,Float64})
        key = SVector(round(pt[1], digits=8), round(pt[2], digits=8), round(pt[3], digits=8))
        idx = get(vertex_map, key, 0)
        if idx == 0
            push!(vertices, pt)
            idx = length(vertices)
            vertex_map[key] = idx
        end
        return idx
    end

    seg_start_vid = Vector{Int}(undef, nseg)
    seg_end_vid   = Vector{Int}(undef, nseg)
    for i in 1:nseg
        seg_start_vid[i] = get_or_add!(sp_arr[i])
        seg_end_vid[i]   = get_or_add!(ep_arr[i])
    end
    empty!(vertex_map)
    GC.gc()

    nv = length(vertices)
    parent_vertex    = zeros(Int, nv)
    children         = [Int[] for _ in 1:nv]
    incoming_segment = zeros(Int, nv)
    for i in 1:nseg
        sv = seg_start_vid[i]
        ev = seg_end_vid[i]
        parent_vertex[ev] = sv
        incoming_segment[ev] = i
        if !(ev in children[sv])
            push!(children[sv], ev)
        end
    end

    root_vertex = seg_start_vid[1]
    return FlowTree(name, vertices, parent_vertex, children, incoming_segment,
                    seg_start_vid, seg_end_vid, seg_diam_cm, seg_label, root_vertex)
end

"""
    load_trees(dict_of_paths) -> Dict{String, FlowTree}

Load multiple trees from a dictionary of `name => csv_path`.
"""
function load_trees(dict_of_paths::Dict{String,String})
    trees = Dict{String, FlowTree}()
    for (name, path) in dict_of_paths
        trees[name] = load_tree(name, path)
        println("  Loaded $(name): $(length(trees[name].segment_start)) segments, $(length(trees[name].vertices)) vertices")
    end
    return trees
end

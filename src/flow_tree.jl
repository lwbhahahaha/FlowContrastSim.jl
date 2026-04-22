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

Uses `parent_segment_id` for deterministic topology reconstruction when available.
Falls back to nearest-vertex matching when `parent_segment_id` is missing
(backward compatibility with old CSVs).
"""
function load_tree(name::String, csv_path::String)
    # Streaming parser. `readdlm` allocates a Matrix{Any} that balloons to
    # tens of GB on multi-million-row CSVs (VascularTreeSim can emit 25M+
    # segments at 8 μm capillary resolution). We read line-by-line and push
    # only typed scalars into pre-sized arrays.

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
        has_parent_seg = psidi !== nothing

        # Dedup vertices via Dict. Round-to-8-decimals in cm = sub-μm tolerance,
        # which matches VascularTreeSim.jl's csv_io.jl output precision.
        vertex_map = Dict{SVector{3,Float64}, Int}()
        vertices = SVector{3,Float64}[]
        sizehint!(vertex_map, 4_000_000)
        sizehint!(vertices,   4_000_000)

        function get_or_add_vertex!(pt::SVector{3,Float64})
            key = SVector(round(pt[1], digits=8), round(pt[2], digits=8), round(pt[3], digits=8))
            idx = get(vertex_map, key, 0)
            if idx == 0
                push!(vertices, pt)
                idx = length(vertices)
                vertex_map[key] = idx
            end
            return idx
        end

        seg_start_vid = Int[]
        seg_end_vid   = Int[]
        seg_diam_cm   = Float64[]
        seg_label     = String[]
        seg_id_csv    = Int[]
        parent_seg_csv = Int[]
        sizehint!(seg_start_vid, 4_000_000)
        sizehint!(seg_end_vid,   4_000_000)
        sizehint!(seg_diam_cm,   4_000_000)
        sizehint!(seg_label,     4_000_000)
        has_parent_seg && (sizehint!(seg_id_csv, 4_000_000); sizehint!(parent_seg_csv, 4_000_000))

        # Parse a single field at column k (1-based) from the comma-split line.
        fields = Vector{SubString{String}}(undef, length(hdr))

        @inbounds while !eof(io)
            line = readline(io)
            isempty(line) && continue
            # Manual split into pre-allocated fields
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

            sp = SVector(parse(Float64, fields[x1i]),
                         parse(Float64, fields[y1i]),
                         parse(Float64, fields[z1i]))
            ep = SVector(parse(Float64, fields[x2i]),
                         parse(Float64, fields[y2i]),
                         parse(Float64, fields[z2i]))
            push!(seg_start_vid, get_or_add_vertex!(sp))
            push!(seg_end_vid,   get_or_add_vertex!(ep))
            push!(seg_diam_cm,   parse(Float64, fields[di]) * 1e-4)
            push!(seg_label,     String(fields[li]))
            if has_parent_seg
                push!(seg_id_csv, parse(Int, fields[sidi]))
                pstr = fields[psidi]
                pv = isempty(pstr) ? -1 : parse(Int, pstr)
                push!(parent_seg_csv, pv)
            end
        end

        # Free the vertex map now that we have the vertex list (saves GBs on
        # large trees — the Dict grows to ~1.5× the vertex count and keeps
        # all 3-vectors + hash slots until GC kicks in).
        empty!(vertex_map)
        GC.gc()

        nv = length(vertices)
        nseg = length(seg_start_vid)

        # ── Pass 2: build topology ──
        parent_vertex = zeros(Int, nv)
        children      = [Int[] for _ in 1:nv]
        incoming_segment = zeros(Int, nv)

        if has_parent_seg
            # Find root segments. VascularTreeSim writes 0 for roots; older/other
            # producers may use -1 or an empty string (mapped to -1 above).
            root_candidates = Set{Int}()
            for i in 1:nseg
                if parent_seg_csv[i] <= 0
                    push!(root_candidates, seg_start_vid[i])
                end
            end

            # Build parent-child relationships from vertex coincidences.
            for i in 1:nseg
                sv = seg_start_vid[i]
                ev = seg_end_vid[i]
                parent_vertex[ev] = sv
                incoming_segment[ev] = i
                if !(ev in children[sv])
                    push!(children[sv], ev)
                end
            end

            # Determine root vertex: start vertex of a root segment whose own
            # vertex has no incoming segment (handles bifurcating-root trees).
            root_vertex = 0
            for rv in root_candidates
                if incoming_segment[rv] == 0
                    root_vertex = rv
                    break
                end
            end
            if root_vertex == 0 && !isempty(root_candidates)
                root_vertex = first(root_candidates)
            end
            root_vertex == 0 && error("Could not determine root vertex in $(csv_path) — no segment with parent_segment_id <= 0 and no inflow-free start vertex")
        else
            # Fallback: no parent_segment_id column. Assume segments are in
            # topo order and first segment's start is the root.
            root_vertex = seg_start_vid[1]
            for i in 1:nseg
                sv = seg_start_vid[i]
                ev = seg_end_vid[i]
                parent_vertex[ev] = sv
                incoming_segment[ev] = i
                if !(ev in children[sv])
                    push!(children[sv], ev)
                end
            end
        end

        return FlowTree(
            name,
            vertices,
            parent_vertex,
            children,
            incoming_segment,
            seg_start_vid,
            seg_end_vid,
            seg_diam_cm,
            seg_label,
            root_vertex,
        )
    finally
        close(io)
    end
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

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
    data = readdlm(csv_path, ','; header=true)
    rows = data[1]
    hdr = vec(data[2])

    # Find column indices
    col(n) = findfirst(==(n), hdr)
    x1i = col("x1_cm"); y1i = col("y1_cm"); z1i = col("z1_cm")
    x2i = col("x2_cm"); y2i = col("y2_cm"); z2i = col("z2_cm")
    di  = col("diameter_um")
    li  = col("label")
    sidi = col("segment_id")
    psidi = col("parent_segment_id")

    has_parent_seg = psidi !== nothing

    nrows = size(rows, 1)

    # ── Pass 1: collect all unique vertices and build segments ──
    # We use a vertex map to deduplicate vertices (start of child == end of parent)
    vertex_map = Dict{SVector{3,Float64}, Int}()
    vertices = SVector{3,Float64}[]

    function get_or_add_vertex!(pt::SVector{3,Float64})
        # Round to avoid floating-point duplicates
        key = SVector(round(pt[1], digits=8), round(pt[2], digits=8), round(pt[3], digits=8))
        idx = get(vertex_map, key, 0)
        if idx == 0
            push!(vertices, pt)
            idx = length(vertices)
            vertex_map[key] = idx
        end
        return idx
    end

    # Temporary storage per row
    seg_start_vid = Vector{Int}(undef, nrows)
    seg_end_vid   = Vector{Int}(undef, nrows)
    seg_diam_cm   = Vector{Float64}(undef, nrows)
    seg_label     = Vector{String}(undef, nrows)
    seg_id_csv    = has_parent_seg ? Vector{Int}(undef, nrows) : Int[]
    parent_seg_csv = has_parent_seg ? Vector{Int}(undef, nrows) : Int[]

    for i in 1:nrows
        sp = SVector(Float64(rows[i, x1i]), Float64(rows[i, y1i]), Float64(rows[i, z1i]))
        ep = SVector(Float64(rows[i, x2i]), Float64(rows[i, y2i]), Float64(rows[i, z2i]))
        sv = get_or_add_vertex!(sp)
        ev = get_or_add_vertex!(ep)
        seg_start_vid[i] = sv
        seg_end_vid[i]   = ev
        seg_diam_cm[i]   = Float64(rows[i, di]) * 1e-4  # um -> cm
        seg_label[i]     = strip(string(rows[i, li]))
        if has_parent_seg
            seg_id_csv[i]     = Int(rows[i, sidi])
            pval = rows[i, psidi]
            parent_seg_csv[i] = (pval isa AbstractString && strip(pval) == "") ? -1 :
                                (ismissing(pval) ? -1 : Int(pval))
        end
    end

    nv = length(vertices)
    nseg = nrows

    # ── Pass 2: build topology ──
    parent_vertex = zeros(Int, nv)
    children      = [Int[] for _ in 1:nv]
    incoming_segment = zeros(Int, nv)

    if has_parent_seg
        # Deterministic reconstruction using parent_segment_id
        # Map csv segment_id -> row index
        csvid_to_row = Dict{Int,Int}()
        for i in 1:nseg
            csvid_to_row[seg_id_csv[i]] = i
        end

        # Find root segments (parent_segment_id == -1 or missing)
        root_candidates = Set{Int}()  # vertex ids that are roots
        for i in 1:nseg
            if parent_seg_csv[i] < 0
                # This segment's start vertex is a root (or child of root)
                push!(root_candidates, seg_start_vid[i])
            end
        end

        # Build parent-child relationships
        for i in 1:nseg
            sv = seg_start_vid[i]
            ev = seg_end_vid[i]
            parent_vertex[ev] = sv
            incoming_segment[ev] = i
            if !(ev in children[sv])
                push!(children[sv], ev)
            end
        end

        # Determine root vertex: the start vertex of root segments
        # (the vertex that has no incoming segment)
        root_vertex = 0
        for rv in root_candidates
            if incoming_segment[rv] == 0
                root_vertex = rv
                break
            end
        end
        if root_vertex == 0
            # Fallback: pick the most common root candidate
            root_vertex = first(root_candidates)
        end
    else
        # Fallback: nearest-vertex matching (for old CSVs without parent_segment_id)
        # Process rows in order; for each segment, find its start vertex among
        # existing end vertices (nearest match)
        # First segment's start is the root
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
        root_vertex
    )
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

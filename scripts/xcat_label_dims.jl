#!/usr/bin/env julia
# Group segments by `label` column and report diameter statistics per label.
# For XCAT segments (label starts with "dias_"), this reveals the anatomical
# diameter profile each tree inherited from NRB. Anomalous thin segments
# (e.g., dias_lad2 at 577 μm nested between 1038/1525 μm neighbors) cause
# downstream flow starvation.

using Printf
using Statistics

const TREE_DIR = length(ARGS) >= 1 ? ARGS[1] :
    "../VascularTreeSim.jl/output"

function label_diam_stats(csv_path::String)
    # Streaming label + diameter read
    io = open(csv_path, "r")
    try
        hdr = split(readline(io), ',')
        li = findfirst(==("label"), hdr)
        di = findfirst(==("diameter_um"), hdr)
        li === nothing && error("no label col")
        di === nothing && error("no diameter_um col")

        label_diams = Dict{String, Vector{Float64}}()
        nline = 0
        while !eof(io)
            line = readline(io); isempty(line) && continue
            nline += 1
            # Only collect XCAT anatomical labels ("dias_*"); grown+sub terminals have generic labels
            # Quick field split
            fields = split(line, ',')
            lbl = String(fields[li])
            # Trim possible trailing CR
            lbl = rstrip(lbl, '\r')
            d = parse(Float64, fields[di])
            push!(get!(label_diams, lbl, Float64[]), d)
        end
        return label_diams
    finally
        close(io)
    end
end

for name in ("LAD", "LCX", "RCA")
    csv = joinpath(TREE_DIR, "$(lowercase(name))_segments.csv")
    isfile(csv) || (println("[$(name)] MISSING"); continue)
    println("=" ^ 80)
    println("[$(name)] CSV: $(csv)")
    t0 = time()
    label_diams = label_diam_stats(csv)
    @printf("  loaded in %.1fs, %d distinct labels\n", time()-t0, length(label_diams))

    # Sort labels: xcat (dias_*) first, then by count desc
    xcat_labels = sort([k for k in keys(label_diams) if startswith(k, "dias_")])
    other_labels = sort([k for k in keys(label_diams) if !startswith(k, "dias_")], by=l -> -length(label_diams[l]))

    println("  --- XCAT anatomical labels (dias_*) ---")
    @printf("  %-25s %10s %12s %12s %12s %12s\n", "label", "count", "min_μm", "median_μm", "max_μm", "mean_μm")
    for lbl in xcat_labels
        d = label_diams[lbl]
        @printf("  %-25s %10d %12.1f %12.1f %12.1f %12.1f\n",
            lbl, length(d), minimum(d), median(d), maximum(d), mean(d))
    end
    println("  --- grown + subdivided labels (sample, top 8 by count) ---")
    @printf("  %-25s %10s %12s %12s %12s %12s\n", "label", "count", "min_μm", "median_μm", "max_μm", "mean_μm")
    for lbl in first(other_labels, 8)
        d = label_diams[lbl]
        @printf("  %-25s %10d %12.1f %12.1f %12.1f %12.1f\n",
            lbl, length(d), minimum(d), median(d), maximum(d), mean(d))
    end
    println()
end

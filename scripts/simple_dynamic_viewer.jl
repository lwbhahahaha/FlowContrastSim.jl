#!/usr/bin/env julia
#
# simple_dynamic_viewer.jl — minimal standalone dynamic-contrast viewer.
# Does NOT use FlowContrastSim's runtime (no full FlowTree, no hemo solver).
#
# Strategy:
#   1. Stream each *_segments.csv, keep only segments with diameter ≥ MIN_DIAM
#      (= 500 µm by default → ~few thousand "main artery" segments / tree).
#   2. Per-segment bolus arrival times come from PHYSICAL hemodynamics:
#      a. If `PEAK_IODINE_DIR/{name}_arrival_time.f32` exists (output of
#         extract_peak_iodine.jl), read it and index by segment_id. This is
#         the real Poiseuille + Pries-viscosity + capillary-bed-R arrival.
#      b. Otherwise fall back to a Murray velocity model: v(s) = v_root · D_s/D_root
#         (volumetric flow scales as D³, area as D², so v scales as D).
#         Better than a global constant velocity, but doesn't account for
#         capillary-bed back-pressure modulation of small-vessel flow.
#   3. Sample a gamma-variate bolus and propagate per segment with the
#      FlowContrastSim dispersion law: disp_factor = sqrt(1 + arrival/t_disp).
#   4. Emit a self-contained Plotly HTML with one Scatter3d trace per tree
#      and a time-slider that flips the marker colors per frame.
#
# Output: OUTPUT_DIR/simple_contrast_viewer.html  (~ a few MB, opens in any browser)
#
# Usage:
#   julia --project=. scripts/simple_dynamic_viewer.jl  [TREE_DIR]  [OUTPUT_DIR]  [MIN_DIAM_UM]  [PEAK_IODINE_DIR]

using Printf

# ── tunables ────────────────────────────────────────────────────────────────
const MIN_DIAM_UM    = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 500.0
const V_ROOT_CM_PER_S = 25.0      # Murray-fallback root velocity (typical proximal LAD/LCX/RCA)
const T_END_S        = 15.0
const N_FRAMES       = 30
const DT_S           = T_END_S / (N_FRAMES - 1)
# Gamma-variate bolus (FlowContrastSim default)
const BOLUS_AMPL  = 5.0
const BOLUS_T0    = 0.5
const BOLUS_TMAX  = 4.0
const BOLUS_ALPHA = 3.0
const T_DISP_S    = 3.0           # bolus dispersion time scale (s); matches FlowContrastSim default
# NRB → phantom-world coord shift (cm), same as elsewhere in the pipeline
const NRB_OFFSET = (2.1443, -9.5553, -20.0068)

# ── args ────────────────────────────────────────────────────────────────────
tree_dir       = length(ARGS) >= 1 ? ARGS[1] : "../VascularTreeSim.jl/output"
out_dir        = length(ARGS) >= 2 ? ARGS[2] : "../phantom_ct_input/contrast_viewer"
peak_iodine_dir = length(ARGS) >= 4 ? ARGS[4] : "../phantom_ct_input/peak_iodine"
isdir(tree_dir) || error("tree_dir not found: $tree_dir")
mkpath(out_dir)

@printf("tree_dir         = %s\n", tree_dir)
@printf("output_dir       = %s\n", out_dir)
@printf("peak_iodine_dir  = %s  %s\n", peak_iodine_dir,
        isdir(peak_iodine_dir) ? "(found)" : "(MISSING → Murray-velocity fallback)")
@printf("min_diam_um      = %.1f\n", MIN_DIAM_UM)
@printf("bolus            = gamma(amp=%.2f mg/mL, t0=%.2fs, tmax=%.2fs, alpha=%.1f)\n",
        BOLUS_AMPL, BOLUS_T0, BOLUS_TMAX, BOLUS_ALPHA)
@printf("t_dispersion_s   = %.2f\n", T_DISP_S)
@printf("n_frames         = %d  (dt=%.3fs, t_end=%.1fs)\n", N_FRAMES, DT_S, T_END_S)
flush(stdout)

# ── gamma-variate bolus + dispersion (same form as FlowContrastSim) ─────────
@inline function gamma_variate(t::Float64)
    t <= BOLUS_T0 && return 0.0
    tp = (t - BOLUS_T0) / (BOLUS_TMAX - BOLUS_T0)
    tp <= 0 && return 0.0
    BOLUS_AMPL * tp^BOLUS_ALPHA * exp(BOLUS_ALPHA * (1.0 - tp))
end
@inline dispersion_div(arrival::Float64) = sqrt(1.0 + arrival / T_DISP_S)

function concentration(t::Float64, arrival::Float64)
    !isfinite(arrival) && return 0.0
    t_shift = t - arrival
    t_shift <= 0 && return 0.0
    disp = dispersion_div(arrival)
    t_eff = BOLUS_T0 + (t_shift - BOLUS_T0) / disp
    t_eff <= BOLUS_T0 && return 0.0
    tp = (t_eff - BOLUS_T0) / (BOLUS_TMAX - BOLUS_T0)
    tp <= 0 && return 0.0
    max(BOLUS_AMPL * tp^BOLUS_ALPHA * exp(BOLUS_ALPHA * (1.0 - tp)) / disp, 0.0)
end

# ── stream a *_segments.csv → kept-segment table ────────────────────────────
struct KeptSeg
    id::Int
    pid::Int
    mx::Float64    # midpoint phantom-world coords (cm)
    my::Float64
    mz::Float64
    len_cm::Float64
    diam_um::Float64
end

function load_main_arteries(csv_path::String)
    ox, oy, oz = NRB_OFFSET
    out = KeptSeg[]
    sizehint!(out, 100_000)
    t0 = time()
    open(csv_path, "r") do io
        readline(io)  # header
        while !eof(io)
            line = readline(io)
            isempty(line) && continue
            # Find first 14 commas (column boundaries)
            c1  = findfirst(',', line);              c2  = findnext(',', line, c1+1)
            c3  = findnext(',', line, c2+1);          c4  = findnext(',', line, c3+1)
            c5  = findnext(',', line, c4+1);          c6  = findnext(',', line, c5+1)
            c7  = findnext(',', line, c6+1);          c8  = findnext(',', line, c7+1)
            c9  = findnext(',', line, c8+1);          c10 = findnext(',', line, c9+1)
            c11 = findnext(',', line, c10+1);         c12 = findnext(',', line, c11+1)
            c13 = findnext(',', line, c12+1);         c14 = findnext(',', line, c13+1)
            # diameter is col 14, cheap parse first to filter
            d_um = parse(Float64, SubString(line, c13+1, c14-1))
            d_um < MIN_DIAM_UM && continue
            sid = parse(Int, SubString(line, c1+1, c2-1))
            pid = parse(Int, SubString(line, c2+1, c3-1))
            x1 = parse(Float64, SubString(line, c3+1, c4-1))
            y1 = parse(Float64, SubString(line, c4+1, c5-1))
            z1 = parse(Float64, SubString(line, c5+1, c6-1))
            x2 = parse(Float64, SubString(line, c6+1, c7-1))
            y2 = parse(Float64, SubString(line, c7+1, c8-1))
            z2 = parse(Float64, SubString(line, c8+1, c9-1))
            len_mm = parse(Float64, SubString(line, c12+1, c13-1))
            push!(out, KeptSeg(
                sid, pid,
                (x1 + x2) / 2.0 + ox,
                (y1 + y2) / 2.0 + oy,
                (z1 + z2) / 2.0 + oz,
                len_mm / 10.0, d_um,
            ))
        end
    end
    @printf("[%s] kept %d segs (diam ≥ %.0f μm)  parse %.1fs\n",
            basename(csv_path), length(out), MIN_DIAM_UM, time()-t0)
    flush(stdout)
    out
end

# ── arrival from f32 file (physical hemo result) ────────────────────────────
# File layout: Float32 vector, length = n_segs_full_tree, indexed 1..n_segs
# by segment_id. We look up each KeptSeg's id directly.
function load_arrival_f32(path::String, segs::Vector{KeptSeg})
    fsz = filesize(path)
    fsz % 4 == 0 || error("$path: size $(fsz) not a multiple of 4 bytes")
    n_file = fsz ÷ 4
    arr_full = Vector{Float32}(undef, n_file)
    open(path, "r") do io
        read!(io, arr_full)
    end
    n = length(segs)
    arrival = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        sid = segs[i].id
        if 1 <= sid <= n_file
            v = arr_full[sid]
            arrival[i] = isfinite(v) ? Float64(v) : Inf
        else
            arrival[i] = Inf
        end
    end
    return arrival
end

# ── Murray-fallback arrival: per-segment velocity scales with diameter ──────
# Murray's law: Q ∝ D³, A ∝ D², so v = Q/A ∝ D. We pick D_root = max diameter
# in the kept set (proxy for the proximal artery) and define v(s) = v_root * D_s/D_root.
# Then transit time τ_s = L_s / v_s, and arrival propagates by topological BFS
# (parent-end → child-start, midpoint convention: arrival = arrival_parent_mid +
# τ_parent/2 + τ_child/2).
function compute_arrival_murray(segs::Vector{KeptSeg})
    n = length(segs)
    n == 0 && return Float64[]
    d_root = maximum(s.diam_um for s in segs)
    @inline function tau_seg(s::KeptSeg)
        v = V_ROOT_CM_PER_S * s.diam_um / d_root
        v <= 0 ? Inf : s.len_cm / v
    end
    id_to_idx = Dict{Int, Int}()
    sizehint!(id_to_idx, n)
    for (i, s) in enumerate(segs)
        id_to_idx[s.id] = i
    end
    children = [Int[] for _ in 1:n]
    roots = Int[]
    for (i, s) in enumerate(segs)
        pi = get(id_to_idx, s.pid, 0)
        if pi > 0
            push!(children[pi], i)
        else
            push!(roots, i)
        end
    end
    arrival = fill(Inf, n)
    queue = Int[]
    sizehint!(queue, n)
    for r in roots
        arrival[r] = tau_seg(segs[r]) / 2.0
        push!(queue, r)
    end
    head = 1
    while head <= length(queue)
        cur = queue[head]; head += 1
        cur_arr = arrival[cur]
        cur_half = tau_seg(segs[cur]) / 2.0
        for ci in children[cur]
            if !isfinite(arrival[ci])
                arrival[ci] = cur_arr + cur_half + tau_seg(segs[ci]) / 2.0
                push!(queue, ci)
            end
        end
    end
    return arrival
end

# Resolve arrival times for a tree: prefer the physical f32 file from
# extract_peak_iodine; fall back to Murray-velocity BFS if it's missing.
function resolve_arrival(segs::Vector{KeptSeg}, tree_name::String, peak_dir::String)
    arr_path = joinpath(peak_dir, lowercase(tree_name) * "_arrival_time.f32")
    if isfile(arr_path)
        @printf("  [%s] arrival source: %s (physical)\n", tree_name, arr_path)
        return load_arrival_f32(arr_path, segs), :physical
    else
        @warn "arrival_time.f32 not found for $tree_name; using Murray-velocity fallback (run extract_peak_iodine.jl first for physical arrival)" arr_path
        return compute_arrival_murray(segs), :murray_fallback
    end
end

# ── discover CSV files (lad/lcx/rca) ────────────────────────────────────────
csvs = sort(filter(f -> endswith(f, "_segments.csv"), readdir(tree_dir; join=true)))
isempty(csvs) && error("No *_segments.csv in $tree_dir")
println("Trees:")
for f in csvs;  println("  $(basename(f))");  end
flush(stdout)

# ── load each tree + compute arrival + per-frame C ──────────────────────────
struct TreeFrame
    name::String
    segs::Vector{KeptSeg}
    arrival::Vector{Float64}
    C::Matrix{Float64}  # n_segs × n_frames
end

times = collect(range(0.0; stop=T_END_S, length=N_FRAMES))
trees = TreeFrame[]
for csv in csvs
    name = uppercase(replace(basename(csv), "_segments.csv" => ""))
    @info "Processing $name"
    segs = load_main_arteries(csv)
    arr, src = resolve_arrival(segs, name, peak_iodine_dir)
    finite_arr = filter(isfinite, arr)
    n_unreach = length(arr) - length(finite_arr)
    if isempty(finite_arr)
        @warn "[$name] no segments with finite arrival — viewer will show all-zero contrast"
    else
        med = sort(finite_arr)[max(1, div(length(finite_arr), 2))]
        @printf("  arrival times (%s): min=%.2fs  median=%.2fs  max=%.2fs  unreachable=%d/%d\n",
                src, minimum(finite_arr), med, maximum(finite_arr),
                n_unreach, length(arr))
    end
    # Per-frame concentration
    C = Matrix{Float64}(undef, length(segs), N_FRAMES)
    for (ti, t) in enumerate(times), (si, a) in enumerate(arr)
        C[si, ti] = concentration(t, a)
    end
    @printf("  C range: %.2f .. %.2f mg/mL  (peak across all time/seg)\n",
            minimum(C), maximum(C))
    push!(trees, TreeFrame(name, segs, arr, C))
    flush(stdout)
end

# ── emit minimal Plotly HTML with time slider ───────────────────────────────
# One Scatter3d trace per tree. Marker mode at segment midpoint, color = C_seg(t).
# Frames update marker.color per trace.

const _BRANCH_COLORS = Dict(
    "LAD" => "#1f77ff", "LCX" => "#e3342f", "RCA" => "#22aa44",
)
branch_color(n) = get(_BRANCH_COLORS, n, "#888888")

# Build base traces — start at first frame
function json_string(s::AbstractString)
    # Minimal escape: just wrap in quotes (assumes no special chars)
    "\"$(replace(s, "\"" => "\\\""))\""
end
function json_floats(v::AbstractVector{<:Real}; digits=4)
    io = IOBuffer()
    print(io, "[")
    for (i, x) in enumerate(v)
        i > 1 && print(io, ",")
        if isfinite(x)
            print(io, round(x; digits=digits))
        else
            print(io, "null")
        end
    end
    print(io, "]")
    String(take!(io))
end

# Per-trace static x/y/z (midpoints) + a representative size (diam → marker)
function build_traces_initial(trees::Vector{TreeFrame})
    io = IOBuffer()
    print(io, "[")
    for (ti, t) in enumerate(trees)
        ti > 1 && print(io, ",")
        xs = [s.mx for s in t.segs]
        ys = [s.my for s in t.segs]
        zs = [s.mz for s in t.segs]
        # marker size: 3..14 pixels scaled by log(diam)
        sizes = [clamp(3 + 11 * log10(max(s.diam_um, MIN_DIAM_UM) / MIN_DIAM_UM) /
                       log10(4000.0 / MIN_DIAM_UM), 3, 14)  for s in t.segs]
        c0 = t.C[:, 1]   # initial frame
        print(io, "{")
        print(io, "\"type\":\"scatter3d\",\"mode\":\"markers\",\"name\":", json_string(t.name), ",")
        print(io, "\"x\":", json_floats(xs; digits=4), ",")
        print(io, "\"y\":", json_floats(ys; digits=4), ",")
        print(io, "\"z\":", json_floats(zs; digits=4), ",")
        print(io, "\"marker\":{")
        print(io, "\"size\":", json_floats(sizes; digits=2), ",")
        print(io, "\"color\":", json_floats(c0; digits=3), ",")
        print(io, "\"colorscale\":[[0,\"#101040\"],[0.05,\"#3050a0\"],[0.4,\"#ffb060\"],[1,\"#ffffe0\"]],")
        print(io, "\"cmin\":0.0,\"cmax\":", round(BOLUS_AMPL; digits=2), ",")
        print(io, "\"colorbar\":", (ti == 1 ? "{\"title\":{\"text\":\"C iodine (mg/mL)\"},\"x\":1.02}" : "null"))
        print(io, "}}")
    end
    print(io, "]")
    String(take!(io))
end

function build_frames_data(trees::Vector{TreeFrame})
    io = IOBuffer()
    print(io, "[")
    for fi in 1:N_FRAMES
        fi > 1 && print(io, ",")
        print(io, "{\"name\":\"f", fi, "\",\"data\":[")
        for (ti, t) in enumerate(trees)
            ti > 1 && print(io, ",")
            ck = t.C[:, fi]
            print(io, "{\"marker\":{\"color\":", json_floats(ck; digits=3), "}}")
        end
        print(io, "]}")
    end
    print(io, "]")
    String(take!(io))
end

function build_slider_steps(times)
    io = IOBuffer()
    print(io, "[")
    for (fi, t) in enumerate(times)
        fi > 1 && print(io, ",")
        print(io, "{\"label\":\"", round(t; digits=2), "s\",")
        print(io, "\"method\":\"animate\",")
        print(io, "\"args\":[[\"f", fi, "\"],{\"mode\":\"immediate\",\"transition\":{\"duration\":0},\"frame\":{\"duration\":0,\"redraw\":true}}]}")
    end
    print(io, "]")
    String(take!(io))
end

# Total segments + axis bounds (for layout)
n_total = sum(length(t.segs) for t in trees)
xmin = minimum(s.mx for t in trees for s in t.segs)
xmax = maximum(s.mx for t in trees for s in t.segs)
ymin = minimum(s.my for t in trees for s in t.segs)
ymax = maximum(s.my for t in trees for s in t.segs)
zmin = minimum(s.mz for t in trees for s in t.segs)
zmax = maximum(s.mz for t in trees for s in t.segs)

@info "Emitting HTML viewer" n_total_segs=n_total n_frames=N_FRAMES

html_path = joinpath(out_dir, "simple_contrast_viewer.html")
open(html_path, "w") do io
    println(io, "<!DOCTYPE html><html lang=\"en\"><head>")
    println(io, "<meta charset=\"utf-8\"><title>Dynamic Contrast — Coronary Trees</title>")
    println(io, "<script src=\"https://cdn.plot.ly/plotly-2.35.2.min.js\" charset=\"utf-8\"></script>")
    println(io, "<style>body{font-family:system-ui,sans-serif;margin:0;padding:8px;background:#101015;color:#eee}h2{margin:4px 0 12px}p{margin:4px 0}#plot{height:90vh}</style>")
    println(io, "</head><body>")
    println(io, "<h2>Dynamic Contrast Transport · LAD + LCX + RCA · diam ≥ $(round(Int, MIN_DIAM_UM)) μm</h2>")
    @printf(io, "<p style=\"color:#aaa;font-size:12px\">Markers at segment midpoints (size ∝ log diameter, color = local iodine in mg/mL). Bolus: gamma-variate amp=%.1f mg/mL t0=%.1fs tmax=%.1fs α=%.0f, dispersion τ=%.1fs. Arrival times from FlowContrastSim hemodynamics if `*_arrival_time.f32` present, else Murray fallback (v_root=%.0f cm/s, v∝D). %d total segments shown across %d frames (Δt=%.2fs).</p>",
            BOLUS_AMPL, BOLUS_T0, BOLUS_TMAX, BOLUS_ALPHA, T_DISP_S, V_ROOT_CM_PER_S,
            n_total, N_FRAMES, DT_S)
    println(io, "<div id=\"plot\"></div>")
    println(io, "<script>")
    println(io, "const data = ", build_traces_initial(trees), ";")
    println(io, "const frames = ", build_frames_data(trees), ";")
    println(io, "const sliderSteps = ", build_slider_steps(times), ";")
    print(io, "const layout = {")
    print(io, "paper_bgcolor:\"#101015\",")
    print(io, "scene:{")
    @printf(io, "xaxis:{title:{text:\"x (cm)\"},range:[%.2f,%.2f],color:\"#ccc\",backgroundcolor:\"#161620\",gridcolor:\"#333\"},",
            xmin - 0.5, xmax + 0.5)
    @printf(io, "yaxis:{title:{text:\"y (cm)\"},range:[%.2f,%.2f],color:\"#ccc\",backgroundcolor:\"#161620\",gridcolor:\"#333\"},",
            ymin - 0.5, ymax + 0.5)
    @printf(io, "zaxis:{title:{text:\"z (cm)\"},range:[%.2f,%.2f],color:\"#ccc\",backgroundcolor:\"#161620\",gridcolor:\"#333\"},",
            zmin - 0.5, zmax + 0.5)
    print(io, "aspectmode:\"data\",bgcolor:\"#101015\"")
    print(io, "},")
    print(io, "margin:{l:0,r:0,b:60,t:30},")
    print(io, "showlegend:true,legend:{font:{color:\"#ccc\"}},")
    print(io, "updatemenus:[{type:\"buttons\",direction:\"left\",x:0.05,y:-0.05,xanchor:\"left\",yanchor:\"top\",pad:{t:10},showactive:false,bgcolor:\"#222\",font:{color:\"#ccc\"},buttons:[{label:\"▶ Play\",method:\"animate\",args:[null,{mode:\"immediate\",fromcurrent:true,frame:{duration:200,redraw:true},transition:{duration:0}}]},{label:\"⏸ Pause\",method:\"animate\",args:[[null],{mode:\"immediate\",frame:{duration:0,redraw:false},transition:{duration:0}}]}]}],")
    print(io, "sliders:[{currentvalue:{prefix:\"t = \",font:{color:\"#fff\"}},pad:{t:30},x:0.12,xanchor:\"left\",len:0.85,bgcolor:\"#222\",activebgcolor:\"#446\",font:{color:\"#ccc\"},steps:sliderSteps}]")
    println(io, "};")
    println(io, "Plotly.newPlot('plot', data, layout, {responsive:true}).then(()=>Plotly.addFrames('plot', frames));")
    println(io, "</script></body></html>")
end

println()
@printf("[viewer] %s (%.1f KB)\n", html_path, filesize(html_path) / 1024)

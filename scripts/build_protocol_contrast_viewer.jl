#!/usr/bin/env julia
#
# build_protocol_contrast_viewer.jl — CTA-like dynamic-contrast HTML viewer for
# the grown LAD/LCX/RCA trees, driven by the SAME injection protocol AIF (Bae
# C_aorta) used elsewhere — not a generic gamma bolus.
#
# Memory-safe: loads + simulates ONE tree at a time, extracts only the
# lightweight viewer subset (top-N segments by diameter + skeleton + per-frame
# concentration), then frees the full ~28 GB FlowTree before the next tree.
#
# Output: OUTPUT_DIR/tree_contrast_viewer.html  (Plotly, time slider, iodine heatmap)
#
# CLI:  julia --project=. scripts/build_protocol_contrast_viewer.jl [TREE_DIR] [OUTPUT_DIR]

using FlowContrastSim
const FCS = FlowContrastSim
using StaticArrays
using LinearAlgebra: norm
using Printf

# ── protocol (saline-chaser x=72, 90 kg) ────────────────────────────────────
const WEIGHT_KG  = 90.0;  const HEIGHT_CM = 173.0
const CONC_MGI   = 370.0
const P1_VOL = 72.0; const P1_RATE = 5.0          # pure iodine
const P2_VOL = 30.0; const P2_RATE = 5.0          # saline chaser (dilution 0)
const KVP = 120.0
const DT = 0.1; const T_END = 35.0                # contrast time grid

# ── hemo + contrast tunables (match extract_peak_iodine_bae) ────────────────
const ROOT_PRESSURE_MMHG     = 100.0
const TERMINAL_PRESSURE_MMHG = 15.0
const HEMATOCRIT             = 0.45
const CAP_BED_R              = 0.12     # hyperemic
const TERRITORY_MASSES_G     = Dict("LAD" => 58.9, "LCX" => 60.9, "RCA" => 63.8)
# ── viewer detail (Standard) ────────────────────────────────────────────────
const MIN_DIAM_UM        = 150.0   # only simulate + show vessels ≥ this
const MAX_SEG_PER_BRANCH = 3000
const TARGET_FRAMES      = 40

const _PALETTE = Dict("LAD"=>"#1f77ff", "LCX"=>"#e3342f", "RCA"=>"#22aa44")
branch_color(n) = get(_PALETTE, n, "#9467bd")

tree_dir = length(ARGS) >= 1 ? ARGS[1] : "../VascularTreeSim.jl/output"
out_dir  = length(ARGS) >= 2 ? ARGS[2] : "../phantom_ct_input/tree_contrast_viewer"
isdir(tree_dir) || error("tree_dir not found: $tree_dir")
mkpath(out_dir)

_interp(times, C, t) = begin
    t <= times[1] && return 0.0
    t >= times[end] && return C[end]
    i = searchsortedlast(times, t); f = (t - times[i]) / (times[i+1] - times[i])
    C[i]*(1-f) + C[i+1]*f
end

# ── 1. Bae AIF for this protocol ────────────────────────────────────────────
@info "Bae PBPK — saline-chaser x=72 (90 kg): phase1 $(P1_VOL)mL@$(P1_RATE) iodine + phase2 $(P2_VOL)mL@$(P2_RATE) saline"
patient  = Patient(weight_kg=WEIGHT_KG, height_cm=HEIGHT_CM)
protocol = TriphasicProtocol(
    weight_kg=WEIGHT_KG, contrast_concentration_mgI_ml=CONC_MGI,
    phase1_volume_ml=P1_VOL, phase1_rate_ml_s=P1_RATE,
    phase2_volume_ml=P2_VOL, phase2_rate_ml_s=P2_RATE,
    phase3_volume_ml=0.0,    phase3_rate_ml_s=5.0,
    phase2_dilution=0.0,     phase3_dilution=0.0)        # phase2 = saline
bae = simulate_central_circulation(patient, protocol; tspan=(0.0, T_END + 5.0), dt_save=0.05)
times_flow = collect(0.0:DT:T_END)
aif_vec = Float64[_interp(bae.times, bae.C_aorta, t) for t in times_flow]
aif_pk = maximum(aif_vec); aif_pt = (argmax(aif_vec) - 1) * DT
@info "AIF (Bae C_aorta): peak $(round(aif_pk,digits=2)) mgI/mL @ t=$(round(aif_pt,digits=1))s; $(length(times_flow)) steps (dt=$DT)"

const STRIDE = max(1, round(Int, length(times_flow) / TARGET_FRAMES))
const FRAME_TIS = collect(1:STRIDE:length(times_flow))
const ALL_TIMES = times_flow[FRAME_TIS]
@info "Frames: $(length(FRAME_TIS)) (stride $STRIDE)"

# ── 2. per-tree lightweight viewer-data extraction (frees tree after) ────────
function extract_viewdata(name, tree, hemo, cr)
    nseg = length(tree.segment_start)
    sparse = !isempty(cr.segment_ids)
    row_of = if sparse
        m = zeros(Int, nseg)
        @inbounds for (i, s) in pairs(cr.segment_ids); 1 <= s <= nseg && (m[s] = i); end
        m
    else
        Int[]
    end
    # Candidate markers = segments that actually carry simulated concentration
    # (sparse set), or any flowing segment in dense mode. Cap UNCONDITIONALLY to
    # the largest-diameter MAX_SEG_PER_BRANCH. (These trees label the bulk
    # "subdivided"/"grown"; we must NOT keep every non-"grown" segment — millions.)
    cands = Int[]
    for s in 1:nseg
        has_data = sparse ? (row_of[s] != 0) : (abs(hemo.segment_flow[s]) > 1e-20)
        has_data && push!(cands, s)
    end
    sort!(cands, by = s -> -tree.segment_diameter_cm[s])
    seg_list = sort(cands[1:min(MAX_SEG_PER_BRANCH, length(cands))])
    nt = length(FRAME_TIS)
    mx=Float64[]; my=Float64[]; mz=Float64[]; diams=Float64[]; lens=Float64[]; ids=Int[]
    conc = [Float64[] for _ in 1:nt]
    for s in seg_list
        a = tree.vertices[tree.segment_start[s]]; b = tree.vertices[tree.segment_end[s]]
        push!(mx, round((a[1]+b[1])/2, digits=3)); push!(my, round((a[2]+b[2])/2, digits=3)); push!(mz, round((a[3]+b[3])/2, digits=3))
        push!(diams, round(1e4*tree.segment_diameter_cm[s], digits=1))
        push!(lens, round(10.0*norm(b-a), digits=3)); push!(ids, s)
        crow = sparse ? row_of[s] : s
        for (fi, ti) in enumerate(FRAME_TIS)
            v = (crow == 0 || ti > size(cr.concentration,2)) ? 0.0 : cr.concentration[crow, ti]
            push!(conc[fi], round(v, digits=2))
        end
    end
    # skeleton lines (subsample whole tree)
    lx=Float64[]; ly=Float64[]; lz=Float64[]
    ls = max(1, Int(ceil(nseg/10000)))
    for s in 1:ls:nseg
        a = tree.vertices[tree.segment_start[s]]; b = tree.vertices[tree.segment_end[s]]
        push!(lx, round(a[1],digits=3)); push!(ly, round(a[2],digits=3)); push!(lz, round(a[3],digits=3))
        push!(lx, round(b[1],digits=3)); push!(ly, round(b[2],digits=3)); push!(lz, round(b[3],digits=3))
        push!(lx, NaN); push!(ly, NaN); push!(lz, NaN)
    end
    sizes = round.([clamp(1.0 + 4.0*sqrt(d/1000.0), 1.5, 8.0) for d in diams], digits=2)
    rf = 0.0
    for c in tree.children[tree.root_vertex]
        sg = tree.incoming_segment[c]; sg != 0 && (rf += hemo.segment_flow[sg])
    end
    cmx = isempty(cr.concentration) ? 0.0 : maximum(cr.concentration)
    return (mx=mx,my=my,mz=mz,lx=lx,ly=ly,lz=lz,diams=diams,lens=lens,ids=ids,
            sizes=sizes,conc=conc,root_flow=rf,cmax=cmx)
end

# ── 3. discover trees + run one at a time ───────────────────────────────────
all_csvs = filter(f -> endswith(lowercase(f), "_segments.csv"), readdir(tree_dir; join=true))
aux = ("domain_points.csv","chambers_points.csv","pericardium_points.csv",
       "great_vessels_points.csv","coronary_arteries_points.csv")
tree_paths = Pair{String,String}[]
for f in sort(all_csvs)
    basename(f) in aux && continue
    m = match(r"^([A-Za-z][A-Za-z0-9]*?)_segments\.csv$", basename(f))
    m === nothing && continue
    push!(tree_paths, uppercase(m.captures[1]) => f)
end
isempty(tree_paths) && error("no *_segments.csv in $tree_dir")
@info "Trees: $(join([n for (n,_) in tree_paths], ", "))"

root_p = ROOT_PRESSURE_MMHG*133.322; term_p = TERMINAL_PRESSURE_MMHG*133.322
branch_data = Dict{String,Any}(); cmax = 0.01
for (name, path) in tree_paths
    t0 = time(); @info "[$name] load_tree …"
    tree = load_tree(name, path)
    @info "[$name] $(length(tree.segment_start)) segs in $(round(time()-t0,digits=1))s"
    hemo = compute_hemodynamics(tree; root_pressure=root_p, terminal_pressure=term_p,
        hematocrit=HEMATOCRIT, capillary_bed_R_per_100g_mmHgmin_ml=CAP_BED_R,
        territory_mass_g=get(TERRITORY_MASSES_G, name, 0.0))
    cr = simulate_contrast(tree, hemo; dt=DT, t_end=T_END, aif=aif_vec, min_diameter_um=MIN_DIAM_UM)
    @info "[$name] contrast peak $(round(maximum(cr.concentration),digits=2)) mgI/mL on $(size(cr.concentration,1)) segs"
    bd = extract_viewdata(name, tree, hemo, cr)
    branch_data[name] = bd; global cmax = max(cmax, bd.cmax)
    @info "[$name] viewer: $(length(bd.mx)) markers  root_flow=$(round(bd.root_flow*60e6,digits=1)) mL/min"
    tree = nothing; hemo = nothing; cr = nothing; GC.gc(); GC.gc()
end

# ── 4. emit HTML (slider + iodine heatmap) ──────────────────────────────────
names = sort(collect(keys(branch_data)))
html = joinpath(out_dir, "tree_contrast_viewer.html")
open(html, "w") do io
    title = "Dynamic Contrast — LAD/LCX/RCA — saline-chaser x=72 (90 kg, 120 kVp), AIF peak $(round(aif_pk,digits=1)) mgI/mL"
    println(io, """<!doctype html><html><head><meta charset="utf-8"><title>$title</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>body{margin:0;font-family:Arial,sans-serif;background:#111;color:#eee}
.controls{padding:8px 14px;display:flex;gap:10px;align-items:center;flex-wrap:wrap;background:#222}
#tr{font-weight:bold;min-width:120px;color:#4af}button{background:#333;color:#eee;border:1px solid #555;padding:4px 12px;cursor:pointer;border-radius:3px}button:hover{background:#444}</style></head><body>
<div style="padding:8px 14px;background:#1a1a1a"><h3 style="margin:0">$title</h3>
<div style="color:#888;font-size:12px">Real Poiseuille+Pries hemodynamics, Taylor-Aris dispersion driven by Bae aortic-root AIF. Color = local iodine (mg/mL).</div></div>
<div class="controls"><button id="pl">▶ Play</button><button id="pa">⏸ Pause</button>
<input id="sl" type="range" min="0" max="$(length(ALL_TIMES)-1)" value="0" step="1" style="width:420px">
<span id="tr">t = 0.00 s</span>
<span style="color:#777;font-size:12px">Root flow: $(join(["$(n)=$(round(branch_data[n].root_flow*60e6,digits=1))" for n in names], ", ")) mL/min</span></div>
<div id="plot" style="width:100vw;height:86vh"></div><script>""")
    println(io, "const allTimes=$(round.(ALL_TIMES,digits=2));")
    println(io, "const traces=[];")
    for (bi,name) in enumerate(names)
        bd = branch_data[name]; col = branch_color(name)
        println(io, "traces.push({type:'scatter3d',mode:'lines',name:'$name',x:$(bd.lx),y:$(bd.ly),z:$(bd.lz),line:{color:'$col',width:1},opacity:0.22,hoverinfo:'skip'});")
        hover = ["$name#$(bd.ids[i]) d=$(bd.diams[i])um" for i in eachindex(bd.ids)]
        println(io, "traces.push({type:'scatter3d',mode:'markers',name:'$name contrast',x:$(bd.mx),y:$(bd.my),z:$(bd.mz),marker:{size:$(bd.sizes),color:$(bd.conc[1]),colorscale:[[0,'#222'],[0.01,'#444'],[0.1,'#2266aa'],[0.4,'#44aaff'],[0.7,'#ffaa22'],[1.0,'#ff2222']],cmin:0,cmax:$(round(cmax,digits=2)),showscale:$(bi==1 ? "true" : "false"),colorbar:{title:'mgI/mL',tickfont:{color:'#aaa'}}},text:$(hover),hoverinfo:'text'});")
        safe = replace(lowercase(name), r"[^a-z0-9]"=>"_")
        println(io, "const $(safe)F=[")
        for fi in 1:length(ALL_TIMES); println(io, fi<length(ALL_TIMES) ? "$(bd.conc[fi])," : "$(bd.conc[fi])"); end
        println(io, "];")
    end
    println(io, "const mIdx=[$(join([2*i-1 for i in 1:length(names)], ","))];")
    println(io, "const aF=[$(join([replace(lowercase(n),r"[^a-z0-9]"=>"_")*"F" for n in names], ","))];")
    println(io, """
function upd(fi){for(let b=0;b<mIdx.length;b++){Plotly.restyle('plot',{'marker.color':[aF[b][fi]]},[mIdx[b]]);}
document.getElementById('tr').textContent='t = '+allTimes[fi].toFixed(2)+' s';document.getElementById('sl').value=fi;}
let pg=false,iv=null;document.getElementById('pl').onclick=()=>{if(pg)return;pg=true;iv=setInterval(()=>{let fi=(Number(document.getElementById('sl').value)+1)%allTimes.length;upd(fi);},120);};
document.getElementById('pa').onclick=()=>{pg=false;if(iv){clearInterval(iv);iv=null;}};
document.getElementById('sl').oninput=e=>upd(Number(e.target.value));
Plotly.newPlot('plot',traces,{scene:{xaxis:{title:'X (cm)',color:'#aaa',gridcolor:'#333'},yaxis:{title:'Y (cm)',color:'#aaa',gridcolor:'#333'},zaxis:{title:'Z (cm)',color:'#aaa',gridcolor:'#333'},aspectmode:'data',bgcolor:'#111'},paper_bgcolor:'#111',margin:{l:0,r:0,b:0,t:0},legend:{font:{color:'#aaa'}}},{displaylogo:false,responsive:true});""")
    println(io, "</script></body></html>")
end
@printf("[viewer] %s  (%.1f MB, %d frames, cmax=%.2f mgI/mL)\n", html, filesize(html)/1e6, length(ALL_TIMES), cmax)

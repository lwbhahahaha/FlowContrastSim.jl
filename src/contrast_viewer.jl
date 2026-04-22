# ── Interactive contrast viewer: time slider, segment coloring ──
# Works with any number of trees with any names (auto-color palette)

# Auto-assign colors from a distinguishable palette
const _VIEWER_PALETTE = [
    "#1f77ff",  # blue
    "#e3342f",  # red
    "#22aa44",  # green
    "#ff7f0e",  # orange
    "#9467bd",  # purple
    "#17becf",  # cyan
    "#d62728",  # dark red
    "#bcbd22",  # olive
    "#e377c2",  # pink
    "#8c564b",  # brown
]

function _branch_color(idx::Int)
    return _VIEWER_PALETTE[mod1(idx, length(_VIEWER_PALETTE))]
end

function build_contrast_viewer(path::AbstractString,
                               trees::Dict{String, FlowTree},
                               hemo_results::Dict{String, HemodynamicsResult},
                               contrast_results::Dict{String, ContrastResult};
                               time_stride::Int=1,
                               title::String="Dynamic Contrast Transport",
                               max_segments_per_branch::Int=8000,
                               branch_colors::Union{Nothing, Dict{String,String}}=nothing)
    # Auto-assign colors if not provided
    branch_names = sort(collect(keys(trees)))
    colors = if branch_colors !== nothing
        branch_colors
    else
        Dict(name => _branch_color(i) for (i, name) in enumerate(branch_names))
    end

    # Common time grid (from first tree)
    first_cr = contrast_results[branch_names[1]]
    all_times = first_cr.times[1:time_stride:end]
    nt = length(all_times)

    # Global max concentration
    cmax = 0.0
    for name in branch_names
        cr = contrast_results[name]
        cmax = max(cmax, maximum(cr.concentration))
    end
    cmax = max(cmax, 0.01)

    # Build per-branch data (only include segments with flow or every Nth for skeleton)
    branch_data = Dict{String, NamedTuple}()
    for name in branch_names
        tree = trees[name]
        cr = contrast_results[name]
        hemo = hemo_results[name]
        nseg = length(tree.segment_start)

        # Map tree segment index → row in cr.concentration. Sparse mode
        # (non-empty segment_ids) means concentration[row_of_seg[s], ti] — and
        # segments absent from segment_ids have no time-resolved data (treated
        # as zero when building conc_frames).
        sparse_mode = !isempty(cr.segment_ids)
        row_of_seg = if sparse_mode
            m = zeros(Int, nseg)
            @inbounds for (i, s) in pairs(cr.segment_ids)
                s in eachindex(m) && (m[s] = i)
            end
            m
        else
            Int[]
        end

        # Select which segments to include:
        # 1. All non-grown segments
        # 2. All grown segments with flow > threshold (and in sparse mode, in the
        #    simulated subset so we actually have concentration data for them)
        included = Set{Int}()
        for s in 1:nseg
            if tree.segment_label[s] != "grown"
                push!(included, s)
            elseif abs(hemo.segment_flow[s]) > 1e-20 &&
                   (!sparse_mode || row_of_seg[s] != 0)
                push!(included, s)
            end
        end

        # If too many, subsample the grown flowing segments by diameter
        if length(included) > max_segments_per_branch
            primary_segs = [s for s in included if tree.segment_label[s] != "grown"]
            grown_segs = [s for s in included if tree.segment_label[s] == "grown"]
            sort!(grown_segs, by=s -> -tree.segment_diameter_cm[s])
            remaining_budget = max_segments_per_branch - length(primary_segs)
            included = Set(primary_segs)
            for s in grown_segs[1:min(remaining_budget, length(grown_segs))]
                push!(included, s)
            end
        end

        # Add skeleton-only segments (every Nth without flow) for visual completeness
        no_flow_segs = [s for s in 1:nseg if !(s in included) && tree.segment_label[s] == "grown"]
        skeleton_stride = isempty(no_flow_segs) ? 1 : max(1, Int(ceil(length(no_flow_segs) / min(2000, length(no_flow_segs)))))
        for i in 1:skeleton_stride:length(no_flow_segs)
            push!(included, no_flow_segs[i])
        end

        seg_list = sort(collect(included))

        mx = Float64[]; my = Float64[]; mz = Float64[]
        diams = Float64[]; lens = Float64[]; seg_ids = Int[]
        conc_frames = [Float64[] for _ in 1:nt]

        for s in seg_list
            a = tree.vertices[tree.segment_start[s]]
            b = tree.vertices[tree.segment_end[s]]
            push!(mx, round((a[1]+b[1])/2, digits=3))
            push!(my, round((a[2]+b[2])/2, digits=3))
            push!(mz, round((a[3]+b[3])/2, digits=3))
            push!(diams, round(1e4 * tree.segment_diameter_cm[s], digits=1))
            push!(lens, round(10.0 * norm(b-a), digits=3))
            push!(seg_ids, s)
            # Row in cr.concentration for this segment; 0 ⇒ simulation skipped it
            crow = sparse_mode ? row_of_seg[s] : s
            for (fi, ti_src) in enumerate(1:time_stride:length(first_cr.times))
                ti_src <= size(cr.concentration, 2) || continue
                fi <= nt || continue
                val = (crow == 0) ? 0.0 : cr.concentration[crow, ti_src]
                push!(conc_frames[fi], round(val, digits=2))
            end
        end

        # Skeleton lines (all segments for visual structure, subsampled)
        lx = Float64[]; ly = Float64[]; lz = Float64[]
        line_stride = max(1, nseg / 10000 |> x -> Int(ceil(x)))
        for s in 1:line_stride:nseg
            a = tree.vertices[tree.segment_start[s]]
            b = tree.vertices[tree.segment_end[s]]
            push!(lx, round(a[1],digits=3)); push!(ly, round(a[2],digits=3)); push!(lz, round(a[3],digits=3))
            push!(lx, round(b[1],digits=3)); push!(ly, round(b[2],digits=3)); push!(lz, round(b[3],digits=3))
            push!(lx, NaN); push!(ly, NaN); push!(lz, NaN)
        end

        sizes = [clamp(1.0 + 4.0 * sqrt(d / 1000.0), 1.5, 8.0) for d in diams]

        branch_data[name] = (
            mx=mx, my=my, mz=mz,
            lx=lx, ly=ly, lz=lz,
            diams=diams, lengths=lens, seg_ids=seg_ids,
            sizes=round.(sizes, digits=2),
            conc_frames=conc_frames,
        )
    end

    total_points = sum(length(bd.mx) for (_, bd) in branch_data)
    println("  Contrast viewer: $(total_points) marker points, $(nt) time frames")

    # ── Generate HTML ──
    open(path, "w") do io
        println(io, """<!doctype html><html><head><meta charset="utf-8"><title>$(title)</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
body{margin:0;font-family:Arial,sans-serif;background:#111}
.controls{padding:8px 14px;display:flex;gap:10px;align-items:center;flex-wrap:wrap;background:#222;color:#eee}
#time-readout{font-weight:bold;min-width:120px;color:#4af}
.info{padding:6px 14px;color:#aaa;font-size:13px;background:#1a1a1a}
button{background:#333;color:#eee;border:1px solid #555;padding:4px 12px;cursor:pointer;border-radius:3px}
button:hover{background:#444}
</style></head><body>
<div style="padding:8px 14px;background:#1a1a1a"><h2 style="margin:0 0 4px;color:#eee">$(title)</h2>
<div style="color:#888;font-size:13px">Poiseuille hemodynamics with Pries viscosity model</div></div>
<div class="controls">
<button id="play-btn">Play</button>
<button id="pause-btn">Pause</button>
<input id="time-slider" type="range" min="0" max="$(nt-1)" value="0" step="1" style="width:400px">
<span id="time-readout">t = 0.00 s</span>
<span style="color:#666;font-size:12px">Root flow: $(join([begin
    tr = trees[name]; hm = hemo_results[name]; rf = 0.0
    for c in tr.children[tr.root_vertex]
        seg = tr.incoming_segment[c]
        seg != 0 && (rf += hm.segment_flow[seg])
    end
    "$(name)=$(round(rf*60e6, digits=1))"
end for name in branch_names], ", ")) mL/min</span>
</div>
<div id="plot" style="width:100vw;height:85vh"></div>
<script>""")

        println(io, "const allTimes = $(round.(all_times, digits=3));")
        println(io, "const cmax = $(round(cmax, digits=3));")

        # Build Plotly traces
        println(io, "const traces = [];")

        for (bi, name) in enumerate(branch_names)
            bd = branch_data[name]
            color = get(colors, name, _branch_color(bi))

            # Skeleton lines
            println(io, "traces.push({type:'scatter3d',mode:'lines',name:'$(name)',x:$(bd.lx),y:$(bd.ly),z:$(bd.lz),line:{color:'$(color)',width:1},opacity:0.25,hoverinfo:'skip'});")

            # Concentration markers
            conc0 = bd.conc_frames[1]
            hover = ["$(name)#$(bd.seg_ids[i]) d=$(bd.diams[i])um L=$(bd.lengths[i])mm" for i in eachindex(bd.seg_ids)]

            println(io, "traces.push({type:'scatter3d',mode:'markers',name:'$(name) contrast',x:$(bd.mx),y:$(bd.my),z:$(bd.mz),marker:{size:$(bd.sizes),color:$(conc0),colorscale:[[0,'#222'],[0.01,'#444'],[0.1,'#2266aa'],[0.4,'#44aaff'],[0.7,'#ffaa22'],[1.0,'#ff2222']],cmin:0,cmax:$(round(cmax,digits=3)),showscale:$(bi==1 ? "true" : "false"),colorbar:{title:'mg/mL',tickfont:{color:'#aaa'},titlefont:{color:'#aaa'}}},text:$(hover),hoverinfo:'text+name'});")

            # Emit frames
            safe_name = replace(lowercase(name), r"[^a-z0-9]" => "_")
            println(io, "const $(safe_name)F=[")
            for fi in 1:nt
                fi < nt ? println(io, "$(bd.conc_frames[fi]),") : println(io, "$(bd.conc_frames[fi])")
            end
            println(io, "];")
        end

        println(io, "const mIdx=[$(join([2*i-1 for i in 1:length(branch_names)], ","))];")
        frame_names = join([replace(lowercase(n), r"[^a-z0-9]" => "_") * "F" for n in branch_names], ",")
        println(io, "const aF=[$(frame_names)];")

        println(io, """
function updateFrame(fi){
  for(let bi=0;bi<mIdx.length;bi++){
    Plotly.restyle('plot',{'marker.color':[aF[bi][fi]]},[mIdx[bi]]);
  }
  document.getElementById('time-readout').textContent='t = '+allTimes[fi].toFixed(2)+' s';
  document.getElementById('time-slider').value=fi;
}
let playing=false,iv=null;
document.getElementById('play-btn').onclick=()=>{
  if(playing)return;playing=true;
  iv=setInterval(()=>{let fi=(Number(document.getElementById('time-slider').value)+1)%allTimes.length;updateFrame(fi);},80);
};
document.getElementById('pause-btn').onclick=()=>{playing=false;if(iv){clearInterval(iv);iv=null;}};
document.getElementById('time-slider').oninput=(e)=>{updateFrame(Number(e.target.value));};
""")

        println(io, "Plotly.newPlot('plot',traces,{scene:{xaxis:{title:'X (cm)',color:'#aaa',gridcolor:'#333'},yaxis:{title:'Y (cm)',color:'#aaa',gridcolor:'#333'},zaxis:{title:'Z (cm)',color:'#aaa',gridcolor:'#333'},aspectmode:'data',bgcolor:'#111'},paper_bgcolor:'#111',plot_bgcolor:'#111',margin:{l:0,r:0,b:0,t:0},legend:{font:{color:'#aaa'}}},{displaylogo:false,responsive:true});")
        println(io, "</script></body></html>")
    end
    return path
end

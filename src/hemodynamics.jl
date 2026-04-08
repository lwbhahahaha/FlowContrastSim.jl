# ── Hemodynamics: Poiseuille flow, resistance, pressure ──
# Works on FlowTree from flow_tree.jl
# Includes Fahraes-Lindqvist effect for non-Newtonian blood in microvessels

const BLOOD_VISCOSITY_PA_S = 0.0035          # 3.5 cP (bulk blood at 45% Hct)
const PLASMA_VISCOSITY_PA_S = 0.0012         # 1.2 cP (plasma only)
const DEFAULT_ROOT_PRESSURE_PA     = 13332.0 # 100 mmHg
const DEFAULT_TERMINAL_PRESSURE_PA = 1999.8  # 15 mmHg
const DEFAULT_DISCHARGE_HEMATOCRIT = 0.45    # 45% systemic hematocrit

struct HemodynamicsResult
    segment_resistance::Vector{Float64}    # Pa*s/m^3
    segment_flow::Vector{Float64}          # m^3/s
    pressure_proximal::Vector{Float64}     # Pa
    pressure_distal::Vector{Float64}       # Pa
    segment_volume_m3::Vector{Float64}     # m^3
    transit_time_s::Vector{Float64}        # seconds
end

# ── Pries et al. (1992, 1994) in-vivo apparent viscosity ──
# Reference: Pries AR, Neuhaus D, Gaehtgens P.
#   "Blood viscosity in tube flow: dependence on diameter and hematocrit"
#   Am J Physiol. 1992;263(6):H1770-8.
# Updated: Pries AR, Secomb TW. Microvascular blood viscosity in vivo
#   and the endothelial surface layer. Am J Physiol. 2005;289:H2657-64.
#
# eta_vitro(D, Hd) models the Fahraes-Lindqvist effect:
#   - For D >> 300 um: eta -> bulk blood viscosity (~3.5 cP at Hd=0.45)
#   - For D ~ 10-300 um: eta decreases (cell-free plasma layer near wall)
#   - For D -> ~2.7 um (RBC diameter): eta rises sharply (RBC squeeze)
#   - D < 2.7 um: flow is impossible
#
# Returns relative viscosity eta_rel = eta_apparent / eta_plasma
function pries_viscosity_relative(diameter_um::Float64;
                                  hematocrit::Float64=DEFAULT_DISCHARGE_HEMATOCRIT)
    D = diameter_um
    Hd = hematocrit

    # Below ~2.7 um, RBCs cannot pass (diameter of RBC ~ 6-8 um but deformable)
    D <= 2.7 && return 1e6  # effectively infinite viscosity

    # Pries 1992 Eq. 1-3 (in-vitro):
    # eta_0.45(D) = 220 * exp(-1.3*D) + 3.2 - 2.44*exp(-0.06*D^0.645)
    eta_045 = 220.0 * exp(-1.3 * D) + 3.2 - 2.44 * exp(-0.06 * D^0.645)

    # C(D) parameter for hematocrit dependence
    # C = (0.8 + exp(-0.075*D)) * (-1 + 1/(1+10^(-11)*D^12)) + 1/(1+10^(-11)*D^12)
    C = (0.8 + exp(-0.075 * D)) * (-1.0 + 1.0 / (1.0 + 1e-11 * D^12)) +
        1.0 / (1.0 + 1e-11 * D^12)

    # eta_rel(D, Hd) = 1 + (eta_0.45 - 1) * ((1-Hd)^C - 1) / ((1-0.45)^C - 1)
    if abs(C) < 1e-10
        eta_rel = eta_045  # degenerate case
    else
        eta_rel = 1.0 + (eta_045 - 1.0) * ((1.0 - Hd)^C - 1.0) / ((1.0 - 0.45)^C - 1.0)
    end

    # In-vivo correction (Pries & Secomb 2005): endothelial surface layer (ESL)
    # The glycocalyx layer (~1.1 um) reduces effective lumen diameter.
    # Only applies for D > 10 um; in smaller vessels the ESL is compressed
    # by passing RBCs and has minimal additional effect.
    if 10.0 < D < 150.0
        w_esl = 1.1  # um, endothelial surface layer thickness
        D_eff = D - 2.0 * w_esl
        # Resistance scales as (D/D_eff)^4 due to Poiseuille
        eta_rel *= (D / D_eff)^4
    end

    return max(eta_rel, 1.0)
end

# Apparent viscosity in Pa*s
function apparent_viscosity(diameter_um::Float64;
                            hematocrit::Float64=DEFAULT_DISCHARGE_HEMATOCRIT)
    return PLASMA_VISCOSITY_PA_S * pries_viscosity_relative(diameter_um; hematocrit=hematocrit)
end

# ── Poiseuille resistance per segment (with diameter-dependent viscosity) ──
function _compute_resistance(tree::FlowTree; hematocrit::Float64=DEFAULT_DISCHARGE_HEMATOCRIT)
    nseg = length(tree.segment_start)
    R = Vector{Float64}(undef, nseg)
    for s in 1:nseg
        a = tree.vertices[tree.segment_start[s]]
        b = tree.vertices[tree.segment_end[s]]
        len_m = norm(b - a) * 0.01       # cm -> m
        rad_m = tree.segment_diameter_cm[s] * 0.5 * 0.01  # cm -> m
        diam_um = tree.segment_diameter_cm[s] * 1e4         # cm -> um
        if rad_m <= 0 || len_m <= 0
            R[s] = 1e30
        else
            mu = apparent_viscosity(diam_um; hematocrit=hematocrit)
            R[s] = 8.0 * mu * len_m / (pi * rad_m^4)
        end
    end
    return R
end

# ── Topological ordering (parent before children, BFS) ──
function _topo_order(tree::FlowTree)
    order = Int[]
    queue = Int[]
    for c in tree.children[tree.root_vertex]
        seg = tree.incoming_segment[c]
        seg != 0 && push!(queue, seg)
    end
    visited = Set{Int}()
    head = 1
    while head <= length(queue)
        s = queue[head]
        head += 1
        s in visited && continue
        push!(visited, s)
        push!(order, s)
        end_v = tree.segment_end[s]
        for c in tree.children[end_v]
            cseg = tree.incoming_segment[c]
            cseg != 0 && push!(queue, cseg)
        end
    end
    return order
end

# ── Identify terminal (leaf) segments ──
function _find_terminals(tree::FlowTree, order::Vector{Int})
    terminals = Set{Int}()
    for s in order
        end_v = tree.segment_end[s]
        has_child = false
        for c in tree.children[end_v]
            cseg = tree.incoming_segment[c]
            if cseg != 0
                has_child = true
                break
            end
        end
        has_child || push!(terminals, s)
    end
    return terminals
end

# ── Subtree resistance (bottom-up) with terminal capillary bed resistance ──
function _subtree_resistance(tree::FlowTree, R::Vector{Float64}, order::Vector{Int};
                             terminal_bed_R::Float64=0.0)
    nseg = length(R)
    sub_R = fill(Inf, nseg)
    terminals = _find_terminals(tree, order)

    # Process in reverse topological order (leaves first)
    for i in length(order):-1:1
        s = order[i]
        end_v = tree.segment_end[s]
        child_segs = Int[]
        for c in tree.children[end_v]
            cseg = tree.incoming_segment[c]
            cseg != 0 && push!(child_segs, cseg)
        end
        if isempty(child_segs)
            # Terminal: segment resistance + capillary bed resistance
            if s in terminals && terminal_bed_R > 0
                sub_R[s] = R[s] + terminal_bed_R
            else
                sub_R[s] = R[s]
            end
        else
            # Parallel children: 1/R_total = sum(1/R_child)
            inv_sum = 0.0
            for cs in child_segs
                isfinite(sub_R[cs]) && (inv_sum += 1.0 / sub_R[cs])
            end
            child_parallel = inv_sum > 0 ? 1.0 / inv_sum : Inf
            sub_R[s] = R[s] + child_parallel
        end
    end
    return sub_R
end

# ── Compute terminal capillary bed resistance for target flow ──
function _calibrate_terminal_resistance(tree::FlowTree, R::Vector{Float64}, order::Vector{Int};
                                         target_flow_m3s::Float64=0.0,
                                         root_pressure::Float64=DEFAULT_ROOT_PRESSURE_PA,
                                         terminal_pressure::Float64=DEFAULT_TERMINAL_PRESSURE_PA)
    target_flow_m3s <= 0 && return 0.0

    dP = root_pressure - terminal_pressure
    terminals = _find_terminals(tree, order)
    n_terminals = length(terminals)
    n_terminals == 0 && return 0.0

    # Binary search for terminal_bed_R that gives target flow
    # First, check flow without terminal resistance
    sub_R0 = _subtree_resistance(tree, R, order; terminal_bed_R=0.0)
    root_segs = Int[]
    for c in tree.children[tree.root_vertex]
        seg = tree.incoming_segment[c]
        seg != 0 && push!(root_segs, seg)
    end
    inv_sum = sum(1.0 / sub_R0[s] for s in root_segs if isfinite(sub_R0[s]))
    flow_no_bed = dP * inv_sum
    flow_no_bed <= target_flow_m3s && return 0.0  # Already below target

    # Binary search
    lo = 0.0
    hi = 1e20  # Very high resistance
    for _ in 1:60
        mid = (lo + hi) / 2
        sub_R_test = _subtree_resistance(tree, R, order; terminal_bed_R=mid)
        inv_s = sum(1.0 / sub_R_test[s] for s in root_segs if isfinite(sub_R_test[s]))
        test_flow = dP * inv_s
        if test_flow > target_flow_m3s
            lo = mid
        else
            hi = mid
        end
    end
    return (lo + hi) / 2
end

# ── Full hemodynamics computation ──
function compute_hemodynamics(tree::FlowTree;
                              root_pressure::Float64=DEFAULT_ROOT_PRESSURE_PA,
                              terminal_pressure::Float64=DEFAULT_TERMINAL_PRESSURE_PA,
                              hematocrit::Float64=DEFAULT_DISCHARGE_HEMATOCRIT,
                              target_flow_ml_min::Float64=0.0)
    nseg = length(tree.segment_start)
    R = _compute_resistance(tree; hematocrit=hematocrit)
    order = _topo_order(tree)

    # Calibrate terminal bed resistance if target flow specified
    terminal_bed_R = 0.0
    if target_flow_ml_min > 0
        target_m3s = target_flow_ml_min / (60.0 * 1e6)
        terminal_bed_R = _calibrate_terminal_resistance(tree, R, order;
            target_flow_m3s=target_m3s,
            root_pressure=root_pressure,
            terminal_pressure=terminal_pressure)
    end

    sub_R = _subtree_resistance(tree, R, order; terminal_bed_R=terminal_bed_R)

    flow = zeros(nseg)
    p_prox = zeros(nseg)
    p_dist = zeros(nseg)

    # Find root segment(s)
    root_segs = Int[]
    for c in tree.children[tree.root_vertex]
        seg = tree.incoming_segment[c]
        seg != 0 && push!(root_segs, seg)
    end

    # Root flow
    if length(root_segs) == 1
        s = root_segs[1]
        flow[s] = (root_pressure - terminal_pressure) / sub_R[s]
        p_prox[s] = root_pressure
        p_dist[s] = root_pressure - flow[s] * R[s]
    else
        for s in root_segs
            flow[s] = isfinite(sub_R[s]) ? (root_pressure - terminal_pressure) / sub_R[s] : 0.0
            p_prox[s] = root_pressure
            p_dist[s] = root_pressure - flow[s] * R[s]
        end
    end

    # Top-down flow distribution (BFS order)
    for s in order
        if !(s in root_segs) && flow[s] == 0.0
            start_v = tree.segment_start[s]
            parent_seg = tree.incoming_segment[start_v]
            if parent_seg != 0
                p_prox[s] = p_dist[parent_seg]
            else
                p_prox[s] = root_pressure
            end
        end

        # Distribute flow to children
        end_v = tree.segment_end[s]
        child_segs = Int[]
        for c in tree.children[end_v]
            cseg = tree.incoming_segment[c]
            cseg != 0 && push!(child_segs, cseg)
        end

        if !isempty(child_segs)
            inv_sum = 0.0
            for cs in child_segs
                isfinite(sub_R[cs]) && (inv_sum += 1.0 / sub_R[cs])
            end
            for cs in child_segs
                if isfinite(sub_R[cs]) && inv_sum > 0
                    flow[cs] = flow[s] * (1.0 / sub_R[cs]) / inv_sum
                end
                p_prox[cs] = p_dist[s]
                p_dist[cs] = p_prox[cs] - flow[cs] * R[cs]
            end
        end
    end

    # Volume and transit time
    vol = Vector{Float64}(undef, nseg)
    tau = Vector{Float64}(undef, nseg)
    for s in 1:nseg
        a = tree.vertices[tree.segment_start[s]]
        b = tree.vertices[tree.segment_end[s]]
        len_m = norm(b - a) * 0.01
        rad_m = tree.segment_diameter_cm[s] * 0.5 * 0.01
        vol[s] = pi * rad_m^2 * len_m
        q = abs(flow[s])
        tau[s] = q > 1e-18 ? vol[s] / q : Inf
    end

    return HemodynamicsResult(R, flow, p_prox, p_dist, vol, tau)
end

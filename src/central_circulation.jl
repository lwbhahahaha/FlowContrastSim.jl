# ── Bae 1998 PBPK model for central-circulation contrast distribution ──
#
# Adapted from Andy's `bae_model_final.jl` Pluto notebook
# (file_from_andy/bae_model_final.jl). Pluto cell markers stripped; logic
# unchanged; small extensions added for our use:
#   - TriphasicProtocol  : 3-phase injection (contrast / contrast / saline)
#   - CentralCirculationResult : packaged C_RV, C_LV, C_aorta, C_pa, C_pv
#   - simulate_central_circulation : top-level convenience wrapper
#   - iodine_to_hu / hu_per_mgI_ml : kVp-dependent HU conversion
#   - bolus_trigger_time, hu_at_time, aif_from_central : SmartPrep-style
#     scan-timing helpers and the AIF bridge to the in-tree Taylor-Aris PDE
#
# Original model: Bae, Heiken, Brink (1998) Radiology 207:647-655.
# Andy's reimplementation reproduces Bae's Figures 4-6 to within ~7-12%
# of his published simulated peaks; mass conservation < 0.01 %.
#
# ~30 compartments / ~50 ODE states. Vessels are well-mixed pools
# (1 ODE: dC/dt = (Q/V) (C_in - C)); Organs add a perfusion-limited
# capillary + ECF pair (2 ODEs with PS = k_ps × Q).
# Topology follows Bae Fig. 3: arm injection → :v_ue_inj → parallel SVC /
# IVC returns → :rh → pulmonary chain (:pulm_in → :lung → :pulm_out) →
# :lh → :aorta_root → systemic organs → venous return.

using OrdinaryDiffEq

# ─────────────────────────────────────────────────────────────────────
# 1. Patient parameters (Bae page 651 scaling formulas)
# ─────────────────────────────────────────────────────────────────────

"""
    Patient(; weight_kg=70.0, height_cm=173.0, sex=:male, cardiac_output_ml_min=NaN)

Patient parameters for the Bae 1998 PBPK model. Blood volume + cardiac
output are computed from `weight_kg` and `height_cm` via Bae's empirical
formulas unless `cardiac_output_ml_min` is set explicitly (use NaN to
auto-compute).
"""
Base.@kwdef struct Patient
    weight_kg::Float64 = 70.0
    height_cm::Float64 = 173.0
    sex::Symbol = :male
    cardiac_output_ml_min::Float64 = NaN
end

function blood_volume_ml(p::Patient)
    W_lb = p.weight_kg * 2.20462
    H_in = p.height_cm / 2.54
    if p.sex == :male
        return 33.164 * H_in^0.725 * W_lb^0.425 - 1229
    else
        return 34.85 * H_in^0.725 * W_lb^0.425 - 1954
    end
end

function cardiac_output_ml_min(p::Patient)
    isnan(p.cardiac_output_ml_min) || return p.cardiac_output_ml_min
    W_lb = p.weight_kg * 2.20462
    H_in = p.height_cm / 2.54
    return 36.36 * H_in^0.725 * W_lb^0.425
end

function patient_scaling_factors(p::Patient)
    # Reference 70-kg adult: BV = 5000 ml, CO = 6500 ml/min
    return (volume = blood_volume_ml(p) / 5000.0,
            flow   = cardiac_output_ml_min(p) / 6500.0)
end

# ─────────────────────────────────────────────────────────────────────
# 2. Compartment topology
# ─────────────────────────────────────────────────────────────────────

Base.@kwdef struct Vessel
    id::Symbol
    volume_ml::Float64
end

Base.@kwdef struct Organ
    id::Symbol
    cap_volume_ml::Float64
    ecf_volume_ml::Float64
end

# Reference (70 kg) compartments. Patient scaling applied at ODE assembly.
const REFERENCE_COMPARTMENTS = [
    # Right heart + pulmonary path
    Vessel(id=:rh,        volume_ml=180.0),
    Vessel(id=:pulm_in,   volume_ml=130.0),
    Organ(id=:lung,       cap_volume_ml=150.0, ecf_volume_ml=144.0),
    Vessel(id=:pulm_out,  volume_ml=160.0),
    Vessel(id=:lh,        volume_ml=180.0),
    # Aortic root — coronary ostium location (AIF for tree PDE)
    Vessel(id=:aorta_root,volume_ml=100.0),
    # Systemic feeders + organs
    Vessel(id=:art_head,     volume_ml=20.0),
    Organ(id=:head,          cap_volume_ml=37.0,  ecf_volume_ml=484.0),
    Vessel(id=:art_ue,       volume_ml=20.0),
    Organ(id=:ue,            cap_volume_ml=12.0,  ecf_volume_ml=2751.0),
    Vessel(id=:art_heart,    volume_ml=20.0),
    Organ(id=:heart_organ,   cap_volume_ml=10.0,  ecf_volume_ml=103.0),
    Vessel(id=:art_bronchus, volume_ml=20.0),
    Organ(id=:bronchus,      cap_volume_ml=5.0,   ecf_volume_ml=144.0),
    Vessel(id=:art_liver,    volume_ml=20.0),  # hepatic artery
    Vessel(id=:art_stspan,   volume_ml=20.0),
    Organ(id=:stspan,        cap_volume_ml=19.0, ecf_volume_ml=112.0),
    Vessel(id=:art_intestine,volume_ml=20.0),
    Organ(id=:intestine,     cap_volume_ml=35.0, ecf_volume_ml=547.0),
    Vessel(id=:art_kidney,   volume_ml=20.0),
    Organ(id=:kidney,        cap_volume_ml=54.0, ecf_volume_ml=89.0),
    Vessel(id=:art_trunk_le, volume_ml=200.0),
    Organ(id=:trunk_le,      cap_volume_ml=57.0, ecf_volume_ml=11002.0),
    Organ(id=:liver,         cap_volume_ml=71.0, ecf_volume_ml=524.0),
    Vessel(id=:portal,       volume_ml=100.0),
    # Venous returns
    Vessel(id=:v_head,       volume_ml=80.0),
    # UE venous return splits into two parallel paths per Bae Fig. 3.
    # Injection enters the 162 ml/min path (antecubital).
    Vessel(id=:v_ue_inj,     volume_ml=40.0),
    Vessel(id=:v_ue_alt,     volume_ml=40.0),
    Vessel(id=:v_heart,      volume_ml=40.0),   # coronary sinus
    Vessel(id=:v_bronchus,   volume_ml=20.0),
    Vessel(id=:v_kidney,     volume_ml=80.0),
    Vessel(id=:v_liver,      volume_ml=100.0),
    Vessel(id=:v_trunk_le_1, volume_ml=2925.0),
    Vessel(id=:v_trunk_le_2, volume_ml=200.0),
    Vessel(id=:ivc_1,        volume_ml=1000.0),
    Vessel(id=:ivc_2,        volume_ml=700.0),
]

# Organ flows (Bae Table 2, ml/min, 70-kg reference)
const ORGAN_FLOWS = Dict{Symbol, Float64}(
    :head => 975.0,
    :ue => 325.0,
    :heart_organ => 260.0,
    :bronchus => 130.0,
    :liver_hepatic_art => 455.0,
    :stspan => 495.0,
    :intestine => 935.0,
    :kidney => 1430.0,
    :trunk_le => 1495.0,
)

# UE venous return split (Bae Fig. 3)
const UE_SPLIT_INJ = 162.0 / 325.0
const UE_SPLIT_ALT = 163.0 / 325.0
const CO_REF       = 6500.0     # ml/min, reference 70-kg adult

# ─────────────────────────────────────────────────────────────────────
# 3. ODE system assembly
# ─────────────────────────────────────────────────────────────────────

struct BaeModelIndex
    vessel_idx::Dict{Symbol, Int}
    organ_cap_idx::Dict{Symbol, Int}
    organ_ecf_idx::Dict{Symbol, Int}
    excreted_idx::Int       # cumulative excreted iodine mass (mg I)
    n_states::Int
end

function build_index(compartments)
    vessel_idx = Dict{Symbol, Int}()
    organ_cap_idx = Dict{Symbol, Int}()
    organ_ecf_idx = Dict{Symbol, Int}()
    i = 0
    for c in compartments
        if c isa Vessel
            i += 1; vessel_idx[c.id] = i
        elseif c isa Organ
            i += 1; organ_cap_idx[c.id] = i
            i += 1; organ_ecf_idx[c.id] = i
        end
    end
    i += 1
    return BaeModelIndex(vessel_idx, organ_cap_idx, organ_ecf_idx, i, i)
end

struct ScaledParams
    V_vessel::Dict{Symbol, Float64}
    V_cap::Dict{Symbol, Float64}
    V_ecf::Dict{Symbol, Float64}
    Q_organ::Dict{Symbol, Float64}     # ml/s (Bae Table 2 × scaling × 1/60)
    CO::Float64                         # ml/s
    k_ps::Float64
    GFR::Float64                        # ml/s (renal filtration)
end

function build_scaled_params(patient::Patient; k_ps::Float64=5.0)
    sf = patient_scaling_factors(patient)
    V_vessel = Dict{Symbol, Float64}()
    V_cap    = Dict{Symbol, Float64}()
    V_ecf    = Dict{Symbol, Float64}()
    for c in REFERENCE_COMPARTMENTS
        if c isa Vessel
            V_vessel[c.id] = c.volume_ml * sf.volume
        else
            V_cap[c.id] = c.cap_volume_ml * sf.volume
            V_ecf[c.id] = c.ecf_volume_ml * sf.volume
        end
    end
    min_to_s = 1.0 / 60.0
    Q_organ = Dict{Symbol, Float64}(
        :head              => ORGAN_FLOWS[:head] * sf.flow * min_to_s,
        :ue                => ORGAN_FLOWS[:ue] * sf.flow * min_to_s,
        :heart_organ       => ORGAN_FLOWS[:heart_organ] * sf.flow * min_to_s,
        :bronchus          => ORGAN_FLOWS[:bronchus] * sf.flow * min_to_s,
        :stspan            => ORGAN_FLOWS[:stspan] * sf.flow * min_to_s,
        :intestine         => ORGAN_FLOWS[:intestine] * sf.flow * min_to_s,
        :kidney            => ORGAN_FLOWS[:kidney] * sf.flow * min_to_s,
        :trunk_le          => ORGAN_FLOWS[:trunk_le] * sf.flow * min_to_s,
        :liver_hepatic_art => ORGAN_FLOWS[:liver_hepatic_art] * sf.flow * min_to_s,
        :liver_total       => (ORGAN_FLOWS[:liver_hepatic_art] +
                                ORGAN_FLOWS[:stspan] +
                                ORGAN_FLOWS[:intestine]) * sf.flow * min_to_s,
        :lung              => CO_REF * sf.flow * min_to_s,
    )
    CO = cardiac_output_ml_min(patient) * min_to_s
    GFR = 0.19 * Q_organ[:kidney] * 0.60   # 19% of renal plasma flow
    return ScaledParams(V_vessel, V_cap, V_ecf, Q_organ, CO, k_ps, GFR)
end

# ── The ODE function ─────────────────────────────────────────────────
# State: u[i] = concentration (mg I/ml) in compartment i, plus one
# cumulative-excreted state.  All flows ml/s; time s; concentration mg/ml.
function bae_ode!(du, u, params, t)
    sp, idx, inj_fn = params

    Cv(id)   = u[idx.vessel_idx[id]]
    Ccap(id) = u[idx.organ_cap_idx[id]]
    Cecf(id) = u[idx.organ_ecf_idx[id]]

    inj_mgI_s = inj_fn(t)   # mg I/s entering antecubital vein (:v_ue_inj)

    # === Right heart === (RA + RV merged; all venous returns join here)
    Q_total = sp.CO
    Q_head_v   = sp.Q_organ[:head]
    Q_ue_inj   = sp.Q_organ[:ue] * UE_SPLIT_INJ
    Q_ue_alt   = sp.Q_organ[:ue] * UE_SPLIT_ALT
    Q_bronchus = sp.Q_organ[:bronchus]
    Q_coronary = sp.Q_organ[:heart_organ]
    Q_ivc_in   = sp.Q_organ[:liver_total] + sp.Q_organ[:kidney] + sp.Q_organ[:trunk_le]
    mass_in_rh = Q_head_v   * Cv(:v_head)     +
                 Q_ue_inj   * Cv(:v_ue_inj)   +
                 Q_ue_alt   * Cv(:v_ue_alt)   +
                 Q_bronchus * Cv(:v_bronchus) +
                 Q_coronary * Cv(:v_heart)    +
                 Q_ivc_in   * Cv(:ivc_2)
    mass_out_rh = Q_total * Cv(:rh)
    du[idx.vessel_idx[:rh]] = (mass_in_rh - mass_out_rh) / sp.V_vessel[:rh]

    # === Pulmonary artery (RV → lung input) ===
    du[idx.vessel_idx[:pulm_in]] = Q_total * (Cv(:rh) - Cv(:pulm_in)) / sp.V_vessel[:pulm_in]

    # === Lung (organ: cap + ECF with PS exchange) ===
    Q_lung = Q_total
    PS_lung = sp.k_ps * Q_lung
    du[idx.organ_cap_idx[:lung]] = (Q_lung * (Cv(:pulm_in) - Ccap(:lung)) -
                                     PS_lung * (Ccap(:lung) - Cecf(:lung))) / sp.V_cap[:lung]
    du[idx.organ_ecf_idx[:lung]] = PS_lung * (Ccap(:lung) - Cecf(:lung)) / sp.V_ecf[:lung]

    # === Pulmonary vein → Left heart → Aorta root ===
    du[idx.vessel_idx[:pulm_out]]  = Q_total * (Ccap(:lung)  - Cv(:pulm_out)) / sp.V_vessel[:pulm_out]
    du[idx.vessel_idx[:lh]]        = Q_total * (Cv(:pulm_out) - Cv(:lh))      / sp.V_vessel[:lh]
    du[idx.vessel_idx[:aorta_root]]= Q_total * (Cv(:lh)       - Cv(:aorta_root)) / sp.V_vessel[:aorta_root]

    # === Arterial feeders (each takes its organ's flow from aorta_root) ===
    for (art_id, organ_flow_id) in (
            (:art_head, :head), (:art_ue, :ue), (:art_heart, :heart_organ),
            (:art_bronchus, :bronchus), (:art_liver, :liver_hepatic_art),
            (:art_stspan, :stspan), (:art_intestine, :intestine),
            (:art_kidney, :kidney), (:art_trunk_le, :trunk_le))
        Q = sp.Q_organ[organ_flow_id]
        du[idx.vessel_idx[art_id]] = Q * (Cv(:aorta_root) - Cv(art_id)) / sp.V_vessel[art_id]
    end

    # === Systemic organs (cap + ECF, no excretion) ===
    function organ_dynamics!(du, organ_id::Symbol, feeder_id::Symbol, Q::Float64)
        PS = sp.k_ps * Q
        Ca = Cv(feeder_id)
        Civ = Ccap(organ_id)
        Cec = Cecf(organ_id)
        du[idx.organ_cap_idx[organ_id]] = (Q * (Ca - Civ) - PS * (Civ - Cec)) / sp.V_cap[organ_id]
        du[idx.organ_ecf_idx[organ_id]] = PS * (Civ - Cec) / sp.V_ecf[organ_id]
    end
    organ_dynamics!(du, :head,        :art_head,        sp.Q_organ[:head])
    organ_dynamics!(du, :ue,          :art_ue,          sp.Q_organ[:ue])
    organ_dynamics!(du, :heart_organ, :art_heart,       sp.Q_organ[:heart_organ])
    organ_dynamics!(du, :bronchus,    :art_bronchus,    sp.Q_organ[:bronchus])
    organ_dynamics!(du, :stspan,      :art_stspan,      sp.Q_organ[:stspan])
    organ_dynamics!(du, :intestine,   :art_intestine,   sp.Q_organ[:intestine])
    organ_dynamics!(du, :trunk_le,    :art_trunk_le,    sp.Q_organ[:trunk_le])

    # === Kidney (with renal excretion via GFR) ===
    Q_k = sp.Q_organ[:kidney]
    PS_k = sp.k_ps * Q_k
    du[idx.organ_cap_idx[:kidney]] = (Q_k * (Cv(:art_kidney) - Ccap(:kidney)) -
                                       PS_k * (Ccap(:kidney) - Cecf(:kidney)) -
                                       sp.GFR * Ccap(:kidney)) / sp.V_cap[:kidney]
    du[idx.organ_ecf_idx[:kidney]] = PS_k * (Ccap(:kidney) - Cecf(:kidney)) / sp.V_ecf[:kidney]
    du[idx.excreted_idx] = sp.GFR * Ccap(:kidney)

    # === Portal vein → Liver (dual supply: hepatic artery + portal) ===
    Q_portal = sp.Q_organ[:stspan] + sp.Q_organ[:intestine]
    portal_in = (sp.Q_organ[:stspan]    * Ccap(:stspan) +
                 sp.Q_organ[:intestine] * Ccap(:intestine)) / Q_portal
    du[idx.vessel_idx[:portal]] = Q_portal * (portal_in - Cv(:portal)) / sp.V_vessel[:portal]

    Q_ha = sp.Q_organ[:liver_hepatic_art]
    Q_total_liver = Q_ha + Q_portal
    liver_in_C = (Q_ha * Cv(:art_liver) + Q_portal * Cv(:portal)) / Q_total_liver
    PS_liv = sp.k_ps * Q_total_liver
    du[idx.organ_cap_idx[:liver]] = (Q_total_liver * (liver_in_C - Ccap(:liver)) -
                                      PS_liv * (Ccap(:liver) - Cecf(:liver))) / sp.V_cap[:liver]
    du[idx.organ_ecf_idx[:liver]] = PS_liv * (Ccap(:liver) - Cecf(:liver)) / sp.V_ecf[:liver]

    # === Venous returns into the right heart (Bae Fig. 3 parallel topology) ===
    du[idx.vessel_idx[:v_head]] = sp.Q_organ[:head] * (Ccap(:head) - Cv(:v_head)) / sp.V_vessel[:v_head]
    # UE: injection enters here on the 162 ml/min path
    du[idx.vessel_idx[:v_ue_inj]] = (Q_ue_inj * Ccap(:ue) + inj_mgI_s -
                                      Q_ue_inj * Cv(:v_ue_inj)) / sp.V_vessel[:v_ue_inj]
    du[idx.vessel_idx[:v_ue_alt]] = Q_ue_alt * (Ccap(:ue) - Cv(:v_ue_alt)) / sp.V_vessel[:v_ue_alt]
    du[idx.vessel_idx[:v_bronchus]] = Q_bronchus * (Ccap(:bronchus) - Cv(:v_bronchus)) / sp.V_vessel[:v_bronchus]
    du[idx.vessel_idx[:v_heart]] = sp.Q_organ[:heart_organ] * (Ccap(:heart_organ) - Cv(:v_heart)) / sp.V_vessel[:v_heart]
    du[idx.vessel_idx[:v_kidney]] = Q_k * (Ccap(:kidney) - Cv(:v_kidney)) / sp.V_vessel[:v_kidney]
    du[idx.vessel_idx[:v_liver]] = Q_total_liver * (Ccap(:liver) - Cv(:v_liver)) / sp.V_vessel[:v_liver]

    # IVC chain
    Q_tle = sp.Q_organ[:trunk_le]
    du[idx.vessel_idx[:v_trunk_le_1]] = Q_tle * (Ccap(:trunk_le) - Cv(:v_trunk_le_1)) / sp.V_vessel[:v_trunk_le_1]
    du[idx.vessel_idx[:v_trunk_le_2]] = Q_tle * (Cv(:v_trunk_le_1) - Cv(:v_trunk_le_2)) / sp.V_vessel[:v_trunk_le_2]
    Q_ivc_total = Q_tle + Q_k + Q_total_liver
    ivc_1_in = (Q_tle * Cv(:v_trunk_le_2) + Q_k * Cv(:v_kidney) +
                Q_total_liver * Cv(:v_liver)) / Q_ivc_total
    du[idx.vessel_idx[:ivc_1]] = Q_ivc_total * (ivc_1_in - Cv(:ivc_1)) / sp.V_vessel[:ivc_1]
    du[idx.vessel_idx[:ivc_2]] = Q_ivc_total * (Cv(:ivc_1) - Cv(:ivc_2)) / sp.V_vessel[:ivc_2]
    return nothing
end

# ─────────────────────────────────────────────────────────────────────
# 4. Injection protocols  (extends AbstractInjectionProtocol from protocol.jl)
# ─────────────────────────────────────────────────────────────────────

"""
    TriphasicProtocol(; phase1_volume_ml, phase1_rate_ml_s, phase2_volume_ml,
                        phase2_rate_ml_s, phase3_volume_ml, phase3_rate_ml_s,
                        contrast_concentration_mgI_ml=370.0,
                        phase2_dilution=1.0, phase3_dilution=0.0)

Three-phase clinical injection (e.g., UCI Cardiac CTP protocol):
- Phase 1: pure contrast, fast loading rate
- Phase 2: contrast at slower rate (or diluted) — keeps AIF plateau
- Phase 3: saline chaser (typically zero dilution)

`phaseN_dilution` is the iodine fraction relative to neat contrast:
1.0 = pure contrast, 0.30 = 30 % contrast + 70 % saline mix, 0.0 = pure
saline. For the UCI Isovue-370 protocol at 70 kg: phase1=49 mL@5 mL/s,
phase2=25 mL@2.5 mL/s (pure contrast, dilution=1.0), phase3=30 mL saline.
"""
Base.@kwdef struct TriphasicProtocol <: AbstractInjectionProtocol
    weight_kg::Float64
    contrast_concentration_mgI_ml::Float64 = 370.0
    phase1_volume_ml::Float64
    phase1_rate_ml_s::Float64
    phase2_volume_ml::Float64
    phase2_rate_ml_s::Float64
    phase3_volume_ml::Float64
    phase3_rate_ml_s::Float64
    phase2_dilution::Float64 = 1.0   # 1.0 = pure contrast in phase 2
    phase3_dilution::Float64 = 0.0   # 0.0 = pure saline (no iodine)
end

function injection_mass_flow(t::Float64, p::TriphasicProtocol)
    t1 = p.phase1_volume_ml / p.phase1_rate_ml_s
    t2 = t1 + p.phase2_volume_ml / p.phase2_rate_ml_s
    t3 = t2 + p.phase3_volume_ml / p.phase3_rate_ml_s
    c0 = p.contrast_concentration_mgI_ml
    if t < 0.0;  return 0.0; end
    if t < t1;   return p.phase1_rate_ml_s * c0; end
    if t < t2;   return p.phase2_rate_ml_s * c0 * p.phase2_dilution; end
    if t < t3;   return p.phase3_rate_ml_s * c0 * p.phase3_dilution; end
    return 0.0
end

function phase_boundary_times(p::TriphasicProtocol)
    ts = Float64[0.0]
    tprev = 0.0
    for (vol, rate) in ((p.phase1_volume_ml, p.phase1_rate_ml_s),
                        (p.phase2_volume_ml, p.phase2_rate_ml_s),
                        (p.phase3_volume_ml, p.phase3_rate_ml_s))
        if vol > 0.0 && rate > 0.0
            tprev += vol / rate
            push!(ts, tprev)
        end
    end
    return ts
end

total_injected_iodine_mg(p::TriphasicProtocol) = (
    p.phase1_volume_ml * p.contrast_concentration_mgI_ml +
    p.phase2_volume_ml * p.contrast_concentration_mgI_ml * p.phase2_dilution +
    p.phase3_volume_ml * p.contrast_concentration_mgI_ml * p.phase3_dilution
)

# ── Bridge: existing UniphaseNoChaser / WithChaser / Biphase variants also
#   feed the Bae model via the same `injection_mass_flow(t, p)` interface.

function injection_mass_flow(t::Float64, p::UniphaseNoChaser)
    duration = p.weight_kg * p.contrast_volume_per_kg / p.injection_rate_ml_s
    t < 0.0 && return 0.0
    t < duration && return p.injection_rate_ml_s * p.contrast_concentration_mgI_ml
    return 0.0
end

function injection_mass_flow(t::Float64, p::UniphaseWithChaser)
    t1 = p.weight_kg * p.contrast_volume_per_kg / p.injection_rate_ml_s
    t2 = t1 + p.weight_kg * p.chaser_volume_per_kg / p.chaser_rate_ml_s
    t < 0.0 && return 0.0
    t < t1 && return p.injection_rate_ml_s * p.contrast_concentration_mgI_ml
    t < t2 && return p.chaser_rate_ml_s * p.contrast_concentration_mgI_ml * p.chaser_dilution
    return 0.0
end

function injection_mass_flow(t::Float64, p::BiphaseNoChaser)
    t1 = p.weight_kg * p.phase1_volume_per_kg / p.phase1_rate_ml_s
    t2 = t1 + p.weight_kg * p.phase2_volume_per_kg / p.phase2_rate_ml_s
    t < 0.0 && return 0.0
    t < t1 && return p.phase1_rate_ml_s * p.contrast_concentration_mgI_ml
    t < t2 && return p.phase2_rate_ml_s * p.contrast_concentration_mgI_ml
    return 0.0
end

function injection_mass_flow(t::Float64, p::BiphaseWithChaser)
    t1 = p.weight_kg * p.phase1_volume_per_kg / p.phase1_rate_ml_s
    t2 = t1 + p.weight_kg * p.phase2_volume_per_kg / p.phase2_rate_ml_s
    t3 = t2 + p.weight_kg * p.chaser_volume_per_kg / p.chaser_rate_ml_s
    t < 0.0 && return 0.0
    t < t1 && return p.phase1_rate_ml_s * p.contrast_concentration_mgI_ml
    t < t2 && return p.phase2_rate_ml_s * p.contrast_concentration_mgI_ml
    t < t3 && return p.chaser_rate_ml_s * p.contrast_concentration_mgI_ml * p.chaser_dilution
    return 0.0
end

# ─────────────────────────────────────────────────────────────────────
# 5. simulate_central_circulation — top-level wrapper
# ─────────────────────────────────────────────────────────────────────

"""
    CentralCirculationResult

Packaged time-series of iodine concentration (mg I/mL) in major
chambers / aorta plus the raw ODE solution for downstream extraction.

Fields:
- `times`          : Vector{Float64} (s)
- `C_aorta`        : aortic root  (= AIF for the in-tree Taylor-Aris PDE)
- `C_RV`           : right heart  (for beam-hardening check)
- `C_LV`           : left heart   (LV blood pool)
- `C_pulm_artery`  : `:pulm_in`
- `C_pulm_vein`    : `:pulm_out`
- `patient`        : input Patient
- `protocol`       : input protocol
- `total_injected_mgI` : reference for mass-balance checks
- `raw_sol`, `raw_idx`, `raw_sp` : ODE solution + index + scaled params
                                    (for advanced extraction of any compartment)
"""
struct CentralCirculationResult
    times::Vector{Float64}
    C_aorta::Vector{Float64}
    C_RV::Vector{Float64}
    C_LV::Vector{Float64}
    C_pulm_artery::Vector{Float64}
    C_pulm_vein::Vector{Float64}
    patient::Patient
    protocol::AbstractInjectionProtocol
    total_injected_mgI::Float64
    raw_sol
    raw_idx::BaeModelIndex
    raw_sp::ScaledParams
end

"""
    simulate_central_circulation(patient, protocol; tspan=(0.0, 120.0), k_ps=5.0, dt_save=0.1)

Solve the Bae 1998 PBPK ODE system for arm-vein injection. Returns a
`CentralCirculationResult` with C_aorta(t), C_RV(t), C_LV(t), etc.

Uses Tsit5 with tight tolerances (rtol=1e-8, atol=1e-11) and explicit
`tstops` at phase boundaries to handle injection-profile discontinuities.

Mass conservation verified: ∫(C_compartment × V) over all compartments
plus cumulative excreted ≈ ∫(injection_mass_flow) within ODE tolerance.
"""
function simulate_central_circulation(patient::Patient,
                                      protocol::AbstractInjectionProtocol;
                                      tspan::Tuple{Float64,Float64}=(0.0, 120.0),
                                      k_ps::Float64=5.0,
                                      dt_save::Float64=0.1)
    idx = build_index(REFERENCE_COMPARTMENTS)
    sp  = build_scaled_params(patient; k_ps=k_ps)
    inj_fn(t) = injection_mass_flow(t, protocol)
    tstops = _protocol_tstops(protocol)
    u0 = zeros(idx.n_states)
    params = (sp, idx, inj_fn)
    prob = ODEProblem(bae_ode!, u0, tspan, params)
    sol = solve(prob, Tsit5(); reltol=1e-8, abstol=1e-11,
                saveat=dt_save, tstops=tstops)
    times = sol.t
    extract(comp) = [u[idx.vessel_idx[comp]] for u in sol.u]
    return CentralCirculationResult(
        times,
        extract(:aorta_root),
        extract(:rh),
        extract(:lh),
        extract(:pulm_in),
        extract(:pulm_out),
        patient,
        protocol,
        _total_iodine(protocol),
        sol, idx, sp,
    )
end

_total_iodine(p::UniphaseNoChaser)   = p.weight_kg * p.contrast_volume_per_kg * p.contrast_concentration_mgI_ml
_total_iodine(p::UniphaseWithChaser) = p.weight_kg * p.contrast_volume_per_kg * p.contrast_concentration_mgI_ml +
                                       p.weight_kg * p.chaser_volume_per_kg   * p.contrast_concentration_mgI_ml * p.chaser_dilution
_total_iodine(p::BiphaseNoChaser)    = p.weight_kg * (p.phase1_volume_per_kg + p.phase2_volume_per_kg) * p.contrast_concentration_mgI_ml
_total_iodine(p::BiphaseWithChaser)  = p.weight_kg * (p.phase1_volume_per_kg + p.phase2_volume_per_kg) * p.contrast_concentration_mgI_ml +
                                       p.weight_kg * p.chaser_volume_per_kg * p.contrast_concentration_mgI_ml * p.chaser_dilution
_total_iodine(p::TriphasicProtocol)  = total_injected_iodine_mg(p)

_protocol_tstops(p::TriphasicProtocol) = phase_boundary_times(p)
_protocol_tstops(p::UniphaseNoChaser) = begin
    duration = p.weight_kg * p.contrast_volume_per_kg / p.injection_rate_ml_s
    [0.0, duration]
end
_protocol_tstops(p::UniphaseWithChaser) = begin
    t1 = p.weight_kg * p.contrast_volume_per_kg / p.injection_rate_ml_s
    t2 = t1 + p.weight_kg * p.chaser_volume_per_kg / p.chaser_rate_ml_s
    [0.0, t1, t2]
end
_protocol_tstops(p::BiphaseNoChaser) = begin
    t1 = p.weight_kg * p.phase1_volume_per_kg / p.phase1_rate_ml_s
    t2 = t1 + p.weight_kg * p.phase2_volume_per_kg / p.phase2_rate_ml_s
    [0.0, t1, t2]
end
_protocol_tstops(p::BiphaseWithChaser) = begin
    t1 = p.weight_kg * p.phase1_volume_per_kg / p.phase1_rate_ml_s
    t2 = t1 + p.weight_kg * p.phase2_volume_per_kg / p.phase2_rate_ml_s
    t3 = t2 + p.weight_kg * p.chaser_volume_per_kg / p.chaser_rate_ml_s
    [0.0, t1, t2, t3]
end

# ─────────────────────────────────────────────────────────────────────
# 6. HU conversion + bolus tracking
# ─────────────────────────────────────────────────────────────────────

"""
    hu_per_mgI_ml(kvp)

Linear enhancement coefficient: HU per (mg I/mL). Calibrated values from
Bae 1998 (25 at 120 kVp) plus iodine k-edge / mass attenuation
literature for other tube voltages.

  kVp  | HU per mgI/mL
  -----+----------------
   80  | 50.0
  100  | 38.0
  120  | 25.0     (Bae's measured 120 kVp value)
  140  | 22.0

Linearly interpolated between table values, clamped at the endpoints.
"""
function hu_per_mgI_ml(kvp::Real)
    if kvp <= 80.0
        return 50.0
    elseif kvp <= 100.0
        return 50.0 + (kvp - 80.0) * (38.0 - 50.0) / 20.0
    elseif kvp <= 120.0
        return 38.0 + (kvp - 100.0) * (25.0 - 38.0) / 20.0
    elseif kvp <= 140.0
        return 25.0 + (kvp - 120.0) * (22.0 - 25.0) / 20.0
    else
        return 22.0
    end
end

iodine_to_hu_delta(c_mgI_ml::Real; kvp::Real=120.0) = hu_per_mgI_ml(kvp) * c_mgI_ml
"""
    iodine_to_hu(c_mgI_ml; baseline_HU=40.0, kvp=120.0)

Convert iodine concentration to absolute HU value, given a baseline
(non-contrast) HU for the tissue / blood pool in question. Typical
baselines: blood ≈ 40 HU, myocardium ≈ 50 HU, fat ≈ -80 HU.
"""
iodine_to_hu(c_mgI_ml::Real; baseline_HU::Real=40.0, kvp::Real=120.0) =
    baseline_HU + iodine_to_hu_delta(c_mgI_ml; kvp=kvp)

"""
    bolus_trigger_time(result; threshold_delta_HU=150.0, kvp=120.0,
                        compartment=:aorta_root)

Simulate SmartPrep-style bolus tracking: return the first time the
selected compartment's HU enhancement (Δ relative to baseline) reaches
`threshold_delta_HU`. Returns `NaN` if the threshold is never met
within the simulation window.

For the UCI protocol the trigger ROI is the descending aorta. We use
`:aorta_root` as the closest available proxy; the descending aorta sees
the same bolus ~1-2 s later (aortic-arch transit) but the SmartPrep
trigger logic uses descending-aorta blood concentration, which is
identical to aortic-root blood concentration once the wave has passed.
The delay matters only at sub-second resolution.
"""
function bolus_trigger_time(result::CentralCirculationResult;
                            threshold_delta_HU::Real=150.0,
                            kvp::Real=120.0,
                            compartment::Symbol=:aorta_root)
    C = compartment === :aorta_root ? result.C_aorta :
        compartment === :rh         ? result.C_RV    :
        compartment === :lh         ? result.C_LV    :
        error("Unsupported compartment: $compartment")
    threshold_c = threshold_delta_HU / hu_per_mgI_ml(kvp)
    idx = findfirst(c -> c >= threshold_c, C)
    return idx === nothing ? NaN : result.times[idx]
end

# ── Linear interpolation on the saved time grid ──
function _interp_grid(times::Vector{Float64}, C::Vector{Float64}, t_query::Real)
    t_query <= times[1]  && return 0.0
    t_query >= times[end] && return C[end]
    i = searchsortedlast(times, t_query)
    frac = (t_query - times[i]) / (times[i+1] - times[i])
    return C[i] * (1.0 - frac) + C[i+1] * frac
end

"""
    chamber_hu_at(result, t, compartment; baseline_HU=40.0, kvp=120.0)

HU value (= baseline + iodine enhancement) of `compartment` at time `t`.
`compartment` ∈ (`:aorta_root`, `:rh`, `:lh`, `:pulm_artery`, `:pulm_vein`).
"""
function chamber_hu_at(result::CentralCirculationResult, t::Real, compartment::Symbol;
                       baseline_HU::Real=40.0, kvp::Real=120.0)
    C = compartment === :aorta_root  ? result.C_aorta :
        compartment === :rh          ? result.C_RV :
        compartment === :lh          ? result.C_LV :
        compartment === :pulm_artery ? result.C_pulm_artery :
        compartment === :pulm_vein   ? result.C_pulm_vein :
        error("Unsupported compartment: $compartment")
    c = _interp_grid(result.times, C, Float64(t))
    return iodine_to_hu(c; baseline_HU=baseline_HU, kvp=kvp)
end

# ─────────────────────────────────────────────────────────────────────
# 7. AIF bridge to the in-tree PDE
# ─────────────────────────────────────────────────────────────────────

"""
    aif_from_central(result; dt=0.1, t_max=nothing) -> (times, aif_mgI_per_ml)

Resample the aortic-root concentration to a uniform grid matching the
in-tree PDE's `dt`/`t_end`. Returns AIF in mgI/mL ready to pass as
`aif=` to `simulate_contrast`.

`t_max` defaults to the last simulated time.
"""
function aif_from_central(result::CentralCirculationResult;
                          dt::Float64=0.1,
                          t_max::Union{Nothing,Float64}=nothing)
    tmax = t_max === nothing ? result.times[end] : t_max
    times = collect(0.0:dt:tmax)
    aif = [_interp_grid(result.times, result.C_aorta, t) for t in times]
    return times, aif
end

# ─────────────────────────────────────────────────────────────────────
# 8. Mass balance helper (diagnostic)
# ─────────────────────────────────────────────────────────────────────

"""
    mass_balance(result)

Returns a NamedTuple with cumulative-injected vs (in-system + excreted)
mass at the final time. Useful for verifying simulation integrity:
`residual_pct` should be < 0.1 % at default solver tolerance.
"""
function mass_balance(result::CentralCirculationResult)
    sp = result.raw_sp
    idx = result.raw_idx
    u_end = result.raw_sol.u[end]
    m = 0.0
    for (id, i) in idx.vessel_idx
        m += sp.V_vessel[id] * u_end[i]
    end
    for (id, i) in idx.organ_cap_idx
        m += sp.V_cap[id] * u_end[i]
    end
    for (id, i) in idx.organ_ecf_idx
        m += sp.V_ecf[id] * u_end[i]
    end
    excreted = u_end[idx.excreted_idx]
    total_injected = result.total_injected_mgI
    residual = total_injected - m - excreted
    return (in_system_mgI = m,
            excreted_mgI  = excreted,
            tracked_mgI   = m + excreted,
            injected_mgI  = total_injected,
            residual_mgI  = residual,
            residual_pct  = 100.0 * residual / total_injected)
end

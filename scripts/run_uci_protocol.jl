#!/usr/bin/env julia
# UCI Cardiac CTP protocol — beam-hardening risk assessment for the
# 49 / 25 / 30 mL triphase injection on a 70 kg patient (GE Revolution,
# 120 kVp, Isovue-370). Uses the Bae 1998 PBPK central-circulation model
# (src/central_circulation.jl) to predict RV and LV iodine concentration
# and HU values at three candidate V2 scan times.
#
# Protocol (per /file_from_andy/UCI Cardiac CTP - ... .pptx
#                + corrected trigger threshold from Wenbo):
#   Phase 1: 49 mL pure Isovue 370 @ 5.0 mL/s   (loading)
#   Phase 2: 25 mL pure Isovue 370 @ 2.5 mL/s   (sustaining — the rate
#                                                 we want to check)
#   Phase 3: 30 mL saline           @ 2.5 mL/s   (chaser)
#   Arm vein, right antecubital.
#
# Scanner: GE Revolution, 120 kVp, descending-aorta SmartPrep ROI.
#   Trigger: baseline + 150 HU enhancement (NOTE: pptx says +80, Wenbo
#            corrected to +150).
#   Scan delay: t_trigger + Δ, Δ ∈ {6.0, 6.1, 6.2} s. Three V2 snapshots.
#
# Beam-hardening criterion: HU_RV at scan time > 300 ⇒ artifact risk.
#
# Outputs:
#   output_uci_protocol/uci_protocol_traces.csv
#       time (s), C_aorta, C_RV, C_LV, C_pulm_art, C_pulm_vein (mgI/mL),
#       HU_aorta, HU_RV, HU_LV (absolute HU at 120 kVp, baseline 40)
#   output_uci_protocol/uci_protocol_summary.txt
#       Trigger time, per-scan-time RV/LV/Aorta HU, BH flag.
#
# Usage:
#   julia --project=. scripts/run_uci_protocol.jl

using FlowContrastSim
using Printf
using Statistics

# ── Inputs ───────────────────────────────────────────────────────────
const WEIGHT_KG   = 70.0
const HEIGHT_CM   = 173.0
const SEX         = :male

const PHASE1_VOL  = 49.0;  const PHASE1_RATE = 5.0
const PHASE2_VOL  = 25.0;  const PHASE2_RATE = 2.5
const PHASE3_VOL  = 30.0;  const PHASE3_RATE = 2.5
const CONC_MGI_ML = 370.0   # Isovue 370

const KVP             = 120.0
const TRIGGER_DELTA   = 150.0   # baseline + 150 HU at descending aorta
const SCAN_DELAYS_S   = (6.0, 6.1, 6.2)
const BASELINE_HU     = 40.0    # blood pool pre-contrast (typical 120 kVp)
const BH_THRESHOLD_HU = 300.0   # RV HU > this ⇒ beam-hardening artifact risk

const TSPAN = (0.0, 120.0)
const DT    = 0.05              # solver save resolution

const OUT_DIR = joinpath(@__DIR__, "..", "output_uci_protocol")
mkpath(OUT_DIR)

# ── Build protocol + run central circulation ─────────────────────────
patient = Patient(weight_kg=WEIGHT_KG, height_cm=HEIGHT_CM, sex=SEX)
protocol = TriphasicProtocol(
    weight_kg                     = WEIGHT_KG,
    contrast_concentration_mgI_ml = CONC_MGI_ML,
    phase1_volume_ml = PHASE1_VOL, phase1_rate_ml_s = PHASE1_RATE,
    phase2_volume_ml = PHASE2_VOL, phase2_rate_ml_s = PHASE2_RATE,
    phase3_volume_ml = PHASE3_VOL, phase3_rate_ml_s = PHASE3_RATE,
    phase2_dilution = 1.0,   # pure contrast in phase 2 (per Wenbo)
    phase3_dilution = 0.0,   # pure saline in phase 3
)

@info "Running Bae PBPK simulation"
@info "  patient: $(WEIGHT_KG) kg male, CO ≈ $(round(FlowContrastSim.cardiac_output_ml_min(patient), digits=0)) mL/min"
@info "  protocol: $(Int(PHASE1_VOL))/$(Int(PHASE2_VOL))/$(Int(PHASE3_VOL)) mL,  $(PHASE1_RATE)/$(PHASE2_RATE)/$(PHASE3_RATE) mL/s"
@info "  total iodine: $(round(total_injected_iodine_mg(protocol), digits=0)) mg"

result = simulate_central_circulation(patient, protocol;
                                       tspan=TSPAN, dt_save=DT)

# Sanity check: mass conservation
mb = mass_balance(result)
@info "  mass balance: injected=$(round(mb.injected_mgI,digits=0)) mg, " *
      "tracked+excreted=$(round(mb.tracked_mgI,digits=0)) mg, " *
      "residual=$(round(mb.residual_pct,digits=4))%"

# ── Bolus tracking + V2 scan times ───────────────────────────────────
t_trigger = bolus_trigger_time(result; threshold_delta_HU=TRIGGER_DELTA, kvp=KVP)
@info "Trigger (DA, +$(TRIGGER_DELTA) HU): t = $(round(t_trigger, digits=2)) s"

scan_results = []
for Δ in SCAN_DELAYS_S
    t_scan = t_trigger + Δ
    hu_rv = chamber_hu_at(result, t_scan, :rh; baseline_HU=BASELINE_HU, kvp=KVP)
    hu_lv = chamber_hu_at(result, t_scan, :lh; baseline_HU=BASELINE_HU, kvp=KVP)
    hu_ao = chamber_hu_at(result, t_scan, :aorta_root; baseline_HU=BASELINE_HU, kvp=KVP)
    bh = hu_rv > BH_THRESHOLD_HU
    push!(scan_results, (delta=Δ, t_scan=t_scan, hu_rv=hu_rv, hu_lv=hu_lv,
                          hu_ao=hu_ao, beam_hardening=bh))
    @printf("  Δ=%.1fs  t_V2=%.2fs  HU_RV=%.0f  HU_LV=%.0f  HU_aorta=%.0f  %s\n",
            Δ, t_scan, hu_rv, hu_lv, hu_ao, bh ? "⚠ BH RISK" : "OK")
end

# ── Peak values across the whole simulation (informational) ───────────
slope = hu_per_mgI_ml(KVP)
peak_RV   = maximum(result.C_RV);    t_peak_RV   = result.times[argmax(result.C_RV)]
peak_LV   = maximum(result.C_LV);    t_peak_LV   = result.times[argmax(result.C_LV)]
peak_AO   = maximum(result.C_aorta); t_peak_AO   = result.times[argmax(result.C_aorta)]
peak_PA   = maximum(result.C_pulm_artery); t_peak_PA = result.times[argmax(result.C_pulm_artery)]
@info "Peak iodine concentrations + HU (informational):"
@printf("  RV:  peak %.2f mgI/mL @ t=%.1fs → HU=%.0f\n",
        peak_RV, t_peak_RV, BASELINE_HU + slope * peak_RV)
@printf("  LV:  peak %.2f mgI/mL @ t=%.1fs → HU=%.0f\n",
        peak_LV, t_peak_LV, BASELINE_HU + slope * peak_LV)
@printf("  AO:  peak %.2f mgI/mL @ t=%.1fs → HU=%.0f\n",
        peak_AO, t_peak_AO, BASELINE_HU + slope * peak_AO)
@printf("  PA:  peak %.2f mgI/mL @ t=%.1fs → HU=%.0f  (right-heart side of pulmonary)\n",
        peak_PA, t_peak_PA, BASELINE_HU + slope * peak_PA)

# ── CSV traces ───────────────────────────────────────────────────────
csv_path = joinpath(OUT_DIR, "uci_protocol_traces.csv")
open(csv_path, "w") do io
    println(io, "t_s,C_aorta_mgI_ml,C_RV_mgI_ml,C_LV_mgI_ml,C_pulm_artery_mgI_ml,C_pulm_vein_mgI_ml," *
                "HU_aorta,HU_RV,HU_LV")
    for i in eachindex(result.times)
        t   = result.times[i]
        c_ao = result.C_aorta[i]; c_rv = result.C_RV[i]; c_lv = result.C_LV[i]
        c_pa = result.C_pulm_artery[i]; c_pv = result.C_pulm_vein[i]
        hu_ao = BASELINE_HU + slope * c_ao
        hu_rv = BASELINE_HU + slope * c_rv
        hu_lv = BASELINE_HU + slope * c_lv
        @printf(io, "%.3f,%.6f,%.6f,%.6f,%.6f,%.6f,%.2f,%.2f,%.2f\n",
                t, c_ao, c_rv, c_lv, c_pa, c_pv, hu_ao, hu_rv, hu_lv)
    end
end
@info "Traces written: $(csv_path)"

# ── Text summary ─────────────────────────────────────────────────────
summary_path = joinpath(OUT_DIR, "uci_protocol_summary.txt")
open(summary_path, "w") do io
    println(io, "UCI Cardiac CTP — beam-hardening risk assessment")
    println(io, "=" ^ 70)
    println(io, "Patient: $(WEIGHT_KG) kg, $(HEIGHT_CM) cm, $(SEX)")
    @printf(io, "  Cardiac output: %.0f mL/min (Bae scaling)\n",
            FlowContrastSim.cardiac_output_ml_min(patient))
    println(io)
    println(io, "Protocol (Isovue $(Int(CONC_MGI_ML))):")
    @printf(io, "  Phase 1: %.1f mL pure contrast @ %.1f mL/s  (duration %.2f s)\n",
            PHASE1_VOL, PHASE1_RATE, PHASE1_VOL/PHASE1_RATE)
    @printf(io, "  Phase 2: %.1f mL pure contrast @ %.1f mL/s  (duration %.2f s)\n",
            PHASE2_VOL, PHASE2_RATE, PHASE2_VOL/PHASE2_RATE)
    @printf(io, "  Phase 3: %.1f mL saline         @ %.1f mL/s  (duration %.2f s)\n",
            PHASE3_VOL, PHASE3_RATE, PHASE3_VOL/PHASE3_RATE)
    @printf(io, "  Total iodine injected: %.0f mg\n",
            total_injected_iodine_mg(protocol))
    println(io)
    println(io, "Scanner (GE Revolution):")
    @printf(io, "  kVp = %.0f,   HU per mg I/mL = %.1f (Bae 1998 fit + literature)\n",
            KVP, slope)
    @printf(io, "  Blood baseline HU = %.0f\n", BASELINE_HU)
    @printf(io, "  Trigger ROI = descending aorta, threshold = baseline + %.0f HU\n",
            TRIGGER_DELTA)
    @printf(io, "  Scan delays after trigger: %s s\n", string(SCAN_DELAYS_S))
    @printf(io, "  Beam-hardening criterion: HU_RV > %.0f at V2 acquisition\n", BH_THRESHOLD_HU)
    println(io)
    println(io, "Mass balance: residual = ", round(mb.residual_pct, digits=4), " %")
    println(io)
    @printf(io, "Trigger time (DA reaches baseline + %.0f HU): t = %.2f s\n",
            TRIGGER_DELTA, t_trigger)
    println(io)
    println(io, "Peak iodine concentrations (whole simulation):")
    @printf(io, "  RV:  %.2f mgI/mL @ t=%.1fs   → HU = %.0f\n", peak_RV, t_peak_RV, BASELINE_HU + slope * peak_RV)
    @printf(io, "  LV:  %.2f mgI/mL @ t=%.1fs   → HU = %.0f\n", peak_LV, t_peak_LV, BASELINE_HU + slope * peak_LV)
    @printf(io, "  AO:  %.2f mgI/mL @ t=%.1fs   → HU = %.0f\n", peak_AO, t_peak_AO, BASELINE_HU + slope * peak_AO)
    @printf(io, "  PA:  %.2f mgI/mL @ t=%.1fs   → HU = %.0f\n", peak_PA, t_peak_PA, BASELINE_HU + slope * peak_PA)
    println(io)
    println(io, "V2 acquisition snapshots:")
    @printf(io, "  %-10s %-10s %-10s %-10s %-10s %s\n",
            "Δ (s)", "t_V2 (s)", "HU_RV", "HU_LV", "HU_aorta", "Beam-Hardening")
    println(io, "  " * "-"^72)
    for r in scan_results
        @printf(io, "  %-10.1f %-10.2f %-10.0f %-10.0f %-10.0f %s\n",
                r.delta, r.t_scan, r.hu_rv, r.hu_lv, r.hu_ao,
                r.beam_hardening ? "⚠ RISK (HU_RV > $(Int(BH_THRESHOLD_HU)))" : "OK")
    end
    println(io)
    println(io, "Model: Bae 1998 PBPK with ~7-16% validation error vs paper")
    println(io, "       (Andy / file_from_andy/bae_model_final.jl).")
    println(io, "       For uniphasic high-rate the model under-predicts by ~16%;")
    println(io, "       worst-case adjustment: multiply HU_RV by 1.19 (= 321.3/270.2)")
    margin_ok_count = count(r -> !r.beam_hardening, scan_results)
    if margin_ok_count == length(scan_results)
        worst_rv = maximum(r.hu_rv for r in scan_results)
        worst_rv_adj = worst_rv * 1.19
        @printf(io, "\n→ Nominal prediction: HU_RV stays at %.0f at all 3 scan times (< %d).\n",
                worst_rv, Int(BH_THRESHOLD_HU))
        @printf(io, "  Worst-case (×1.19) upper bound: HU_RV ≈ %.0f → still %s threshold.\n",
                worst_rv_adj, worst_rv_adj > BH_THRESHOLD_HU ? "EXCEEDS" : "below")
        println(io, "→ Phase 2 @ 2.5 mL/s does NOT trigger beam-hardening on this model.")
    else
        println(io, "\n→ Beam-hardening RISK at one or more scan times. Consider lower phase-2 rate.")
    end
end
@info "Summary written: $(summary_path)"
println("\nDone.")

#!/usr/bin/env julia
# Find phase 2 (rate, scan-delay) combinations that satisfy:
#   1.  RV HU at scan time ∈ [180, 200]      (target window)
#   2.  RV HU at scan time < 300              (no beam-hardening)
#   3.  Aorta/LV HU at scan time ≥ 350       (typical CTA target)
#
# Holds phase 1 (49 mL @ 5 mL/s) and phase 3 (30 mL saline @ 2.5 mL/s) fixed.
# Sweeps phase 2 rate ∈ {0.5..3.0 mL/s} × scan delay after trigger
# ∈ {6.0..14.0 s}.
#
# Output: output_uci_protocol/sweep_rv_target.csv + console table.

using FlowContrastSim
using Printf

const WEIGHT_KG  = 70.0
const PHASE1_VOL = 49.0;  const PHASE1_RATE = 5.0
const PHASE2_VOL = 25.0
const PHASE3_VOL = 30.0;  const PHASE3_RATE = 2.5
const CONC       = 370.0
const KVP        = 120.0
const TRIGGER_HU = 150.0
const BASELINE   = 40.0
const RV_LOW, RV_HIGH = 180.0, 200.0
const BH_THR = 300.0
const AORTA_MIN_CTA = 350.0

const PHASE2_RATES = [0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0]
const SCAN_DELAYS  = [6.0, 8.0, 10.0, 12.0, 14.0]

const OUT_DIR = joinpath(@__DIR__, "..", "output_uci_protocol")
mkpath(OUT_DIR)
csv_path = joinpath(OUT_DIR, "sweep_rv_target.csv")

patient = Patient(weight_kg=WEIGHT_KG, height_cm=173.0)
slope = hu_per_mgI_ml(KVP)

println("UCI Cardiac CTP — RV target 180-200 HU + CTA quality search")
println("Constants: phase1 49 mL@5 mL/s,  phase3 30 mL saline@2.5 mL/s")
println("Sweep: phase2 rate × scan_delay → RV/LV/Aorta HU at scan time\n")

# Header
@printf("%-6s | %-8s | %-7s | %-7s | %-7s | %-7s | %s\n",
        "rate", "Δ_scan", "t_scan", "HU_RV", "HU_LV", "HU_aorta", "flags")
println("-"^85)

results = []

open(csv_path, "w") do io
    println(io, "phase2_rate_ml_s,phase2_duration_s,scan_delay_s,trigger_t_s,t_scan_s,"
                * "HU_RV,HU_LV,HU_aorta,HU_PA,HU_PV,"
                * "RV_in_target,BH_safe,CTA_aorta_ok")
    for rate in PHASE2_RATES
        prot = TriphasicProtocol(
            weight_kg=WEIGHT_KG, contrast_concentration_mgI_ml=CONC,
            phase1_volume_ml=PHASE1_VOL, phase1_rate_ml_s=PHASE1_RATE,
            phase2_volume_ml=PHASE2_VOL, phase2_rate_ml_s=rate,
            phase3_volume_ml=PHASE3_VOL, phase3_rate_ml_s=PHASE3_RATE,
            phase2_dilution=1.0, phase3_dilution=0.0,
        )
        result = simulate_central_circulation(patient, prot; tspan=(0.0, 120.0), dt_save=0.05)
        t_trig = bolus_trigger_time(result; threshold_delta_HU=TRIGGER_HU, kvp=KVP)
        ph2_dur = PHASE2_VOL / rate
        for Δ in SCAN_DELAYS
            t_scan = t_trig + Δ
            t_scan > result.times[end] && continue
            hu_rv = chamber_hu_at(result, t_scan, :rh;          baseline_HU=BASELINE, kvp=KVP)
            hu_lv = chamber_hu_at(result, t_scan, :lh;          baseline_HU=BASELINE, kvp=KVP)
            hu_ao = chamber_hu_at(result, t_scan, :aorta_root;  baseline_HU=BASELINE, kvp=KVP)
            hu_pa = chamber_hu_at(result, t_scan, :pulm_artery; baseline_HU=BASELINE, kvp=KVP)
            hu_pv = chamber_hu_at(result, t_scan, :pulm_vein;   baseline_HU=BASELINE, kvp=KVP)
            rv_ok  = RV_LOW <= hu_rv <= RV_HIGH
            bh_ok  = hu_rv < BH_THR
            cta_ok = hu_ao >= AORTA_MIN_CTA
            flags  = string(rv_ok ? "RV★" : "   ", " ", bh_ok ? "BH✓" : "BH✗", " ",
                            cta_ok ? "CTA✓" : "CTA✗")
            push!(results, (; rate, Δ, t_scan, hu_rv, hu_lv, hu_ao, rv_ok, bh_ok, cta_ok))
            @printf("%-6.2f | %-8.1f | %-7.2f | %-7.0f | %-7.0f | %-7.0f | %s\n",
                    rate, Δ, t_scan, hu_rv, hu_lv, hu_ao, flags)
            @printf(io, "%.3f,%.4f,%.1f,%.4f,%.4f,%.2f,%.2f,%.2f,%.2f,%.2f,%s,%s,%s\n",
                    rate, ph2_dur, Δ, t_trig, t_scan, hu_rv, hu_lv, hu_ao, hu_pa, hu_pv,
                    rv_ok ? "yes" : "no", bh_ok ? "yes" : "no", cta_ok ? "yes" : "no")
        end
    end
end

# Summary: which combinations hit RV target AND CTA quality
println("\n" * "="^85)
println("RV ∈ [180, 200] HU candidates (★ rows):")
for r in results
    if r.rv_ok
        flag = r.cta_ok ? "✓ CTA aorta ≥ $(Int(AORTA_MIN_CTA)) HU" : "✗ aorta only $(round(Int,r.hu_ao)) HU (too low for CTA)"
        @printf("  rate=%.2f mL/s,  Δ_scan=%.1fs  →  RV=%.0f,  LV=%.0f,  aorta=%.0f  %s\n",
                r.rate, r.Δ, r.hu_rv, r.hu_lv, r.hu_ao, flag)
    end
end
println()
println("CSV: $(csv_path)")

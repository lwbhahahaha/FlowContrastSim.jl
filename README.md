# FlowContrastSim.jl

Hemodynamic flow simulation and dynamic contrast transport on vascular trees
produced by [VascularTreeSim.jl](../VascularTreeSim.jl) (or any CSV that
follows the same schema).

Given a tree CSV plus boundary pressures + a per-100 g downstream resistance,
the package computes Poiseuille flow with diameter-dependent non-Newtonian
viscosity (Pries 1992 / Pries & Secomb 2005), distributes flow across the
tree, and (optionally) simulates a contrast bolus propagating through the
network with a per-segment time-of-arrival + dispersion model. Results are
exported to a self-contained 3-D HTML viewer.

The simulator is **physics-only — no calibration**. Root flow is whatever the
tree topology, BCs, viscosity, and the literature-derived capillary-bed
resistance produce; nothing back-solves to hit a target. The package is also
**organ-agnostic** — pass a brain-vessel or skeletal-muscle CSV with matching
configs and the same pipeline runs.

---

## Pipeline

```
<tree>_segments.csv  (at-rest geometry)        <tree>_segments.csv  (max-dilated)
            │                                              │
            ▼                                              ▼
┌─────────────────────────┐                  ┌─────────────────────────┐
│ coronary_baseline.toml  │                  │ coronary_hyperemic.toml │
│   cap_R = 0.15          │                  │   cap_R = 0.12          │
│   (Pries-Secomb at rest)│                  │   (= 0.15 × autoreg     │
│                         │                  │    cap-bed relax 20 %)  │
└─────────────────────────┘                  └─────────────────────────┘
            │                                              │
            ▼                                              ▼
┌─────────────────────────┐                  ┌─────────────────────────┐
│ run_flow_simulation     │                  │ run_flow_simulation     │
│   1. load_tree (parent- │                  │      (same code path)   │
│      id topology)       │                  │                         │
│   2. Poiseuille + Pries │                  │                         │
│   3. R_terminal_bed     │                  │                         │
│      from cap_R × mass  │                  │                         │
│   4. top-down flow      │                  │                         │
│   5. contrast transport │                  │                         │
└─────────────────────────┘                  └─────────────────────────┘
            │                                              │
            ▼                                              ▼
       baseline ~58 (LAD)                        hyperemic ~184 (LAD)
       CFR ≈ 3-4× per tree

   Add `scripts/lad_stenosis_autoreg_sweep.jl` to fold closed-loop
   arteriolar autoregulation on top: arterioles dilate up to 1.6×
   (Wong-Molloi 2008 reserve) to maintain target Q as stenosis grows;
   when reserve is exhausted, Q falls.
```

---

## Key concepts

**load_tree uses `parent_segment_id`, not coordinate dedup.** Dense subdivided
trees (~110 M segments / LCX in ~60 g of myocardium) pack endpoints tightly
enough that 10 nm coordinate rounding can merge unrelated vertices. The older
dedup-by-coordinates loader overwrote `incoming_segment[v]` and orphaned a few
thousand grown segments together with ~23 M subdivided descendants from the
flow BFS — observed as "21% of LCX segs not flowing". The current loader
walks the `parent_segment_id` chain that the grower writes deterministically:
each segment gets its own unique end vertex, its start vertex is its parent's
end vertex. No coordinate coincidence is consulted. The bug is permanently
gone; running it on a brain tree, a muscle tree, or any other CSV that
follows the schema works the same way.

**Pries + Fahraeus-Lindqvist viscosity.** Apparent viscosity tracks Pries 1992
(in-vitro) corrected by Pries & Secomb 2005's endothelial surface layer
(adds an effective glycocalyx in 10–150 μm vessels):
- bulk blood viscosity at large D
- decreases at 10–300 μm (cell-free plasma layer near wall)
- rises sharply below ~10 μm (RBC squeeze) — divergent at 2.7 μm

**Physics-grounded downstream R, no target back-solve.** Each terminal
segment gets a series resistance such that the per-territory total is the
literature value:

```
R_per_terminal = (capillary_bed_R_per_100g_mmHgmin_ml × 100 / territory_mass_g) × n_terminals
```

When all `n_terminals` are combined in parallel back at the root, the
effective downstream R equals the territory's per-100 g cap-bed R. The
hyperemic value 0.15 mmHg·min/mL/100 g comes directly from Pries & Secomb
2005 for myocardium. The same parameter applies to baseline — the state
difference comes from the **tree geometry**, not from changing this knob (see
"strict dilation" below).

`compute_hemodynamics` also still accepts `target_flow_ml_min` to back-solve a
uniform terminal-bed R via binary search. **Don't use that on production runs**
— it violates the no-calibration rule and is kept only for legacy benchmarks.

**Strict dilation per Wong-Molloi 2008.** Two coronary states use two trees:

- **Max-dilated**: the as-grown tree (Murray-optimal diameters). Used for the
  hyperemic config.
- **At-rest**: same tree with arterioles in `[8, 400]` μm scaled by 0.625
  (= 1/1.6, reversing Wong-Molloi 2008's empirical 1.6× max-hyperemia
  multiplier). See `VascularTreeSim.jl/scripts/scale_to_rest.jl`. Conduits
  (> 400 μm, structural rigidity) and capillaries (< 8 μm, no SM) unchanged.

cap_R is `0.15` at rest (Pries-Secomb in-vivo cap-bed + venous return) and
`0.12` under hyperemia (the cap bed itself shows a small autoregulatory
relaxation under max vasodilation, on top of the arteriolar dilation).
CFR (coronary flow reserve) is then a property of the **anatomical scaling
+ a small cap-bed knob**, not a free fit to a flow target — matching how
real myocardium controls perfusion. This also gives the right vessel
diameters to a virtual CT scanner: an at-rest CT sees constricted
arterioles, a hyperemic CT sees max-dilated arterioles, both from the same
underlying anatomy.

**Closed-loop autoregulation** (`scripts/lad_stenosis_autoreg_sweep.jl`)
goes one step further: instead of two fixed geometries, it solves for the
arteriolar dilation factor f ∈ [0.625, 1.0] that maintains a target root
flow as stenosis is applied. Uses an analytical R-model (R_total =
R_other(s) + b/f^4) fit from two no-stenosis hemos, so each later stenosis
level needs only 2-3 hemo calls (vs ~8 for bisection). When the reserve is
exhausted (f hits 1.0 and Q still < target), the simulation flags it and
reports the residual flow.

**Sparse contrast.** The default `contrast_min_diameter_um = 50` skips
capillary-level segments from the `(n_segs × n_timesteps)` concentration
matrix. Without this a 110 M-segment LCX would need ~50 GB just for the
contrast field. Hemodynamics (R + flow + transit time) is still computed for
every segment.

**Real physical arrival times (no log-compression).** Bolus arrival at
each segment is the hemodynamics-derived `vol/flow` cumulative transit
(via topological BFS) — no truncation or rescaling. Previous releases
applied a "log-compress to [0, max_arrival_s]" hack for visualization
which contaminated the peak-time estimate and downstream voxelization;
that hack is removed. `max_arrival_s` is still accepted as a kwarg for
API compatibility but silently ignored.

**Taylor-Aris analytical PDE solution (when an AIF is provided).** With a
protocol-derived AIF on the root, the contrast field on every segment is
the exact path-integrated Green's-function solution of the 1D advection-
dispersion equation:

```
C_s(t)   = (AIF ⊛ G_{μ_s, σ_s})(t)
μ_s     = Σ_{ancestors} τ_a   + τ_s / 2
σ_s²    = Σ_{ancestors} σ_a²  + σ_seg_s² / 2
σ_seg²  = 2 D_eff(R, v) × L / v³
D_eff   = D_mol + (R v)² / (48 D_mol)        # Taylor 1953 / Aris 1956
```

Per-segment mean and variance accumulate together in a single BFS pass.
Mass conservation is exact (∫C_s(t) dt = ∫AIF dt for every reachable
segment) because the Gaussian convolution is unit-area.

What this physics gives you that the empirical `sqrt(1 + a/t_dispersion_s)`
formula didn't (= why we bothered to swap):

- **No tuned dispersion constant.** Dispersion is derived from the
  segment's radius, length, velocity, and the iodine molecular diffusivity
  `D_mol_m2_s` (default 1.5e-9, iohexol/iomeprol in plasma at 37 °C). The
  knob `t_dispersion_s` becomes legacy — only consulted when no AIF is
  supplied. Aligns with the project's "no calibration" rule.
- **Different protocols produce physically correct downstream
  differences.** A sharper injection (high rate, low volume) gives a
  sharper AIF; the new path preserves that sharpness in proximal segments
  and smears it in deeper segments by an amount determined by the actual
  tree geometry. The old `sqrt` formula erased protocol-to-protocol
  differences after a few seconds of transit, regardless of geometry.
- **Per-tree and per-state dispersion emerges automatically.** LAD, LCX,
  RCA have different path-length distributions → their effective AIF →
  tissue transfer functions differ. At-rest vs hyperemic geometries
  (Wong-Molloi scaling) shift arteriolar velocities → Taylor variance
  scales correctly. The previous `t_dispersion_s` was a single global
  number that ignored both effects.
- **Closed-loop validation with perfusion_pipeline becomes meaningful.**
  When `basis_simulator` renders volumes and the perfusion pipeline
  re-measures AIF + tissue curves to recover MBF, the relationship
  between them is now governed by true physics (anatomy + diffusion). With
  an empirical dispersion knob, any recovered "transit dispersion" would
  just be the knob value rather than physiology.
- **Stenosis sweeps make sense.** Stenosis changes velocities by orders
  of magnitude in the affected branch. Taylor variance scales as
  `~v²·L/v³ = L/v`, so dispersion responds correctly to velocity changes.
  `sqrt(1 + a/t_disp)` was nearly velocity-independent and gave the same
  smear regardless of obstruction.

What it does *not* change: root flow (set by hemodynamics, not contrast),
arrival-time means (pure advection, already exact in the previous code),
or peak time at the root segment (≈ AIF peak time, regardless of model).
The differences accumulate downstream.

The legacy gamma-variate path with `contrast_t_dispersion_s` is still the
default when no `[injection_protocol]` is present in the config — for
backward compatibility with old configs and ablation studies.

**`ContrastResult.arrival_s`.** Each `ContrastResult` now carries a
length-`n_seg` `arrival_s::Vector{Float64}` of physical arrival times.
`extract_peak_iodine.jl` writes this as a Float32 `{tree}_arrival_time.f32`
alongside `{tree}_peak_iodine.f32` — consumed by the perfusion pipeline's
`make_basissim_phantoms.jl --use-myo-arrival` and the standalone
`simple_dynamic_viewer.jl`.

**Protocol-driven AIF (mgI/mL).** The legacy `(contrast_t0, contrast_tmax,
contrast_amplitude, contrast_alpha)` knobs are a hand-tuned gamma-variate
with no physical units. The new path in `src/protocol.jl` synthesizes the
AIF at the aorta root from a clinical injection protocol +
patient-physiology parameters:

```
protocol (vol/rate/conc)  →  injection_profile(t)     mgI/s into peripheral vein
                          ⊛  central_transit_kernel    gamma-variate impulse response
                                                       of RV → lung → LV → aorta
                          ÷  cardiac_output_ml_s       → AIF(t) in mgI/mL
```

The central-transit kernel is a *lumped phenomenological* gamma-variate
(peak time + dispersion), not a chamber-by-chamber mechanistic model —
matches the agreed scope: we do not solve PDEs in the chambers, only in
the three coronary trees. Mass conservation is exact: ∫AIF dt =
total_injected_iodine / cardiac_output.

Currently supported protocols (subtypes of `AbstractInjectionProtocol`):

| protocol              | phases                                                       | extra fields beyond uniphase                                                                       |
|-----------------------|--------------------------------------------------------------|----------------------------------------------------------------------------------------------------|
| `UniphaseNoChaser`    | one rectangular pulse                                        | (none)                                                                                              |
| `UniphaseWithChaser`  | contrast pulse → chaser (saline or mixed)                    | `chaser_volume_per_kg`, `chaser_rate_ml_s`, `chaser_dilution`                                       |
| `BiphaseNoChaser`     | two contrast pulses at different rates                       | `phase1_volume_per_kg`, `phase1_rate_ml_s`, `phase2_volume_per_kg`, `phase2_rate_ml_s`              |
| `BiphaseWithChaser`   | two contrast pulses → chaser                                 | both biphase fields + the three chaser fields above                                                 |

All four share `weight_kg` and `contrast_concentration_mgI_ml` (default
370 mgI/mL ≈ Iomeron 370). `chaser_dilution` is the iodine fraction of
the chaser (0 = pure saline, 0.30 = 30 % contrast / 70 % saline mixed
bolus, 1.0 = pure contrast). All defaults are TCGA-typical clinical
values.

Phase concatenation uses **proportional boundary blending** (the sample
at a phase boundary is the area-weighted mix of the two phases' fluxes
over its `dt`-wide interval), so the rectangle-rule integral is exact:
`∫ injection_profile dt = total_injected_iodine` regardless of whether
each phase's duration is an integer multiple of `dt`.

Adding more schemes (e.g. arterial-phase bolus shaping for stenosis
studies) is three small edits: define a new `<: AbstractInjectionProtocol`
struct, list its phases as `(volume_ml, rate_ml_s, concentration_mgI_ml)`
tuples in an `injection_profile(p)` method, and add one `elseif` to
`_parse_injection_protocol`. `protocol_to_aif`, `central_transit_kernel`,
and the Taylor-Aris solver dispatch generically and do not need editing.

When `[injection_protocol]` is present in the TOML (or a protocol is
passed to `run_flow_simulation` via kwarg), AIF is built once and shared
across all three trees, then fed as the root boundary condition to the
**Taylor-Aris analytical PDE solver** described above. The legacy
gamma-variate knobs and the empirical `sqrt(1+a/t_disp)` dispersion
become fallback only — used when the protocol section is absent. AIF
unit is **mgI/mL** consistently through the pipeline, so the contrast
field can flow directly into `basis_simulator` material decomposition
without rescaling.

---

## Installation

Requires Julia ≥ 1.9. Not registered in General; develop from a local clone:

```julia
using Pkg
Pkg.develop(path="/path/to/FlowContrastSim.jl")
Pkg.instantiate()
```

Dependencies are `StaticArrays`, `TOML`, `NearestNeighbors`, plus Julia
stdlibs. No GPU dependency.

---

## Quick start — coronary

Assume you have grown trees in `../VascularTreeSim.jl/output/` (max-dilated)
and `../VascularTreeSim.jl/output_at_rest/` (post-`scale_to_rest.jl`). The
canonical natural-physics report is:

```bash
# baseline (at-rest tree)
julia --project=. scripts/natural_flow_summary.jl \
      ../VascularTreeSim.jl/output_at_rest \
      configs/coronary_baseline.toml

# hyperemic (max-dilated tree)
julia --project=. scripts/natural_flow_summary.jl \
      ../VascularTreeSim.jl/output \
      configs/coronary_hyperemic.toml
```

Each invocation loads three trees (~3 min each) and prints per-tree root
flow + % segments flowing. Full report ~15 min per state.

For the dynamic contrast pipeline (writes a Plotly HTML viewer):

```julia
using FlowContrastSim
cfg = load_flow_config("configs/coronary_hyperemic.toml")
result = run_flow_simulation("../VascularTreeSim.jl/output", cfg;
                             output_dir="out_hyperemic")
# result.trees, result.hemo_results, result.contrast_results, result.viewer_path
```

`scripts/lad_stenosis_sweep.jl` shrinks every `dias_lad1` (proximal LAD)
segment by 0–90% in 10 % increments and writes a CSV / PNG of the Gould-style
flow vs stenosis curve.

### Dynamic contrast with a clinical injection protocol

`configs/coronary_hyperemic_uniphase.toml` runs the same hyperemic
physics with the AIF now derived from a single-phase, no-chaser
injection protocol (see the "Protocol-driven AIF" key concept). The
`[injection_protocol]` and `[patient_physiology]` blocks are the
user-facing knobs; everything else stays the same:

```julia
using FlowContrastSim
cfg = load_flow_config("configs/coronary_hyperemic_uniphase.toml")
result = run_flow_simulation("../VascularTreeSim.jl/output", cfg;
                             output_dir="out_uniphase")
# result.aif, result.protocol, result.physiology are also returned.
```

Three ways to set the protocol variables per run:

```julia
# 1) Edit the TOML, re-load.

# 2) Pass overrides as kwargs to run_flow_simulation — kwargs win over TOML.
run_flow_simulation("../VascularTreeSim.jl/output", cfg;
    injection_protocol = UniphaseNoChaser(
        weight_kg                     = 90.0,
        contrast_concentration_mgI_ml = 320.0,
        contrast_volume_per_kg        = 0.6,
        injection_rate_ml_s           = 6.0),
    patient_physiology = PatientPhysiology(
        cardiac_output_ml_s          = 100.0,    # e.g. younger patient
        central_transit_delay_s      = 10.0,
        central_transit_dispersion_s = 2.5))

# 3) Build a modified config once, run multiple times against it.
cfg2 = with_protocol(cfg; protocol=UniphaseNoChaser(weight_kg=90.0))
run_flow_simulation("../VascularTreeSim.jl/output", cfg2)
```

The same patterns work with the biphase/chaser protocols:

```julia
# Biphase fast-loading + slow-sustaining, with chaser
run_flow_simulation("../VascularTreeSim.jl/output", cfg;
    injection_protocol = BiphaseWithChaser(
        weight_kg              = 80.0,
        phase1_volume_per_kg   = 0.4,   # fast loading
        phase1_rate_ml_s       = 6.0,
        phase2_volume_per_kg   = 0.4,   # slow sustain
        phase2_rate_ml_s       = 3.0,
        chaser_volume_per_kg   = 0.5,
        chaser_rate_ml_s       = 3.0,
        chaser_dilution        = 0.30))  # 30:70 mix
```

### Peak-only / SVP snapshot mode

For Single-Volume Perfusion (SVP) downstream — where the simulator only
needs to feed `basis_simulator` a single iodine snapshot at the V2
acquisition time, not the full time series — pass `peak_time_s` to
either `simulate_contrast` or `run_flow_simulation`:

```julia
# Use AIF peak + 2 s as the V2 acquisition time (typical bolus-tracking trigger delay)
_, aif = protocol_to_aif(cfg.injection_protocol, cfg.patient_physiology;
                          dt=cfg.dt, t_max=cfg.t_end)
v2_t = (argmax(aif) - 1) * cfg.dt + 2.0    # ≈ 17.8 s for default uniphase

result = run_flow_simulation("../VascularTreeSim.jl/output", cfg; peak_time_s=v2_t)

# result.contrast_results["LAD"].concentration is now (n_sim_segs × 1),
# .times = [v2_t]. Memory shrinks ~N_t× compared to the full time series.
```

This avoids allocating the `(n_sim_segs × n_t)` per-tree matrix and runs
exactly one Gaussian convolution sample per segment instead of computing
the full `AIF ⊛ G_σ` on the time grid.

### Just the AIF (no flow simulation)

```julia
times, aif = protocol_to_aif(cfg.injection_protocol, cfg.patient_physiology;
                              dt=cfg.dt, t_max=cfg.t_end)
# aif :: Vector{Float64}, mgI/mL, on the same time grid as the simulation
```

`examples/uniphase_protocol_aif.jl` is the end-to-end demo and writes
`output/uniphase_aif.csv` for quick plotting. The diagnostic script
`scripts/run_8um_uniphase.jl` is the full 8 μm phantom-grown tree run
(loads 100 M+ segments per tree, runs hemodynamics + Taylor-Aris PDE in
sparse mode at ≥50 μm + peak-only snapshot, writes
`output_8um_uniphase/summary.txt`).

---

## Config TOML schema

```toml
root_pressure_mmhg    = 100.0
terminal_pressure_mmhg = 15.0
discharge_hematocrit  = 0.45

# physics-grounded downstream R (Pries-Secomb in-vivo for myocardium):
capillary_bed_R_per_100g_mmHgmin_ml = 0.15   # 0.12 in coronary_hyperemic.toml

# Legacy contrast bolus shape (gamma-variate). Used only when no
# [injection_protocol] block is present; otherwise the protocol-derived AIF
# overrides these fields. See "Protocol-driven AIF" above. For real patient
# AIF measurements, see ../perfusion_pipeline/scripts/step0_prepare_aif.py
contrast_amplitude = 5.0   # mg I / mL peak
contrast_t0   = 0.5        # s onset
contrast_tmax = 4.0        # s peak
contrast_alpha = 3.0       # shape
contrast_t_dispersion_s = 3.0   # bolus dispersion timescale: disp = sqrt(1 + a/t_disp)
dt    = 0.1                # s — simulation timestep
t_end = 20.0               # s
max_arrival_s = 15.0       # LEGACY — accepted but ignored after the log-compress fix
contrast_min_diameter_um = 50.0   # skip contrast on segs below this diameter

# Iodine molecular diffusivity (m²/s) used in the Taylor-Aris D_eff term
# (only consulted on the AIF / PDE path). Default is iohexol/iomeprol in
# plasma at 37 °C. Sensitivity is logarithmic (the Taylor term and the
# molecular term enter with opposite D_mol-dependence), so this rarely
# needs changing.
D_mol_m2_s = 1.5e-9

# per-tree territory mass for cap_R parallel→series conversion
[territory_masses_g]
LAD = 58.9
LCX = 60.9
RCA = 63.8

# clinical reference values for the SUMMARY printout (display only — they
# do NOT enter the physics)
[target_flows_ml_min]
LAD = 50.0
LCX = 30.0
RCA = 50.0

# ── Optional: protocol-driven AIF (mgI/mL).
# When present, replaces the legacy gamma-variate root input. See
# `configs/coronary_hyperemic_uniphase.toml` for a full example.
[injection_protocol]
type                          = "UniphaseNoChaser"  # currently: UniphaseNoChaser
weight_kg                     = 70.0
contrast_concentration_mgI_ml = 370.0   # default 370 (Iomeron 370 / Omnipaque 350-equivalents)
contrast_volume_per_kg        = 0.5     # default 0.5 mL/kg
injection_rate_ml_s           = 5.0     # default 5 mL/s

# ── Optional: patient physiology (only consulted when [injection_protocol]
# is set). Defaults are healthy resting adult; literature ranges:
#   cardiac_output_ml_s          50–117    (3–7 L/min)
#   central_transit_delay_s      8–16
#   central_transit_dispersion_s 2–5
[patient_physiology]
cardiac_output_ml_s          = 83.0    # = 5 L/min
central_transit_delay_s      = 12.0
central_transit_dispersion_s = 3.0
```

`configs/coronary_baseline.toml` uses `cap_R = 0.15` (at-rest, Pries-Secomb)
paired with `output_at_rest/`; `configs/coronary_hyperemic.toml` uses
`cap_R = 0.12` (autoregulatory cap-bed relaxation) paired with `output/`
(max-dilated). `configs/coronary_hyperemic_uniphase.toml` is the same
hyperemic physics but uses the protocol-driven AIF with a uniphase
injection. `configs/coronary_hyperemic_biphase_chaser.toml` shows the
biphase + saline-chaser CTCA-style protocol. To target a brain tree,
swap the masses, swap the tree dir, leave the physics knobs alone.

---

## CSV contract

Input columns (produced by `VascularTreeSim.jl`):

| column | unit | role |
|---|---|---|
| `segment_id` | int | unique within tree |
| `parent_segment_id` | int | **authoritative topology** (0 = root, > 0 = id of segment whose end vertex is this segment's start vertex) |
| `x1,y1,z1, x2,y2,z2` | cm | endpoint coords (used only for length computation) |
| `length_mm`, `diameter_um` | — | scalar per-segment fields |
| `label` | — | `dias_*` for XCAT, `grown` / `subdivided` for generated branches |

The loader rebuilds topology from `parent_segment_id` only. Older trees
without that column fall through to a legacy coordinate-dedup path; emit a
warning if you see "detached roots".

---

## Output

`run_flow_simulation` returns `(trees, hemo_results, contrast_results,
viewer_path)`. `hemo_results[name]` is a `HemodynamicsResult` with
per-segment vectors:

| field | unit | what |
|---|---|---|
| `segment_resistance` | Pa·s/m³ | Poiseuille × Pries viscosity |
| `segment_flow` | m³/s | top-down distributed flow |
| `pressure_proximal`, `pressure_distal` | Pa | per-segment pressures |
| `segment_volume_m3` | m³ | π r² L |
| `transit_time_s` | s | volume / |flow| |

`contrast_results[name]` is a `ContrastResult` with the following fields
(all sparse when `contrast_min_diameter_um > 0`):

| field | dimensions | what |
|---|---|---|
| `times` | `n_timesteps` (or `1` in peak-only mode) | s |
| `concentration`, `outlet_concentration` | `n_sim_segs × n_timesteps` | mgI/mL |
| `segment_ids` | `n_sim_segs` | indices into the full tree; empty in dense mode (`row i ↔ segment i`) |
| `arrival_s` | `n_seg` full tree | physical bolus arrival at midpoint (s); `Inf` for unreachable segments |
| `arrival_variance_s2` | `n_seg` full tree | Taylor-Aris cumulative variance σ² (s²) at midpoint; populated on the AIF / PDE path, `Inf` on the legacy gamma path |

`arrival_variance_s2` is the per-segment dispersion handle the perfusion
pipeline can use to flag distal segments where the Gaussian convolution
broadens the bolus by more than `~ k × dt` (i.e., where the σ-window
covers many time samples). Useful when comparing simulated vs measured
tissue curve widths in the downstream MBF deconvolution.

The Plotly HTML viewer (`viewer_path`) animates this with a time slider +
play/pause + per-segment hover.

---

## Diagnostic scripts

| script | purpose |
|---|---|
| `natural_flow_summary.jl` | minimal hemo-only report for a directory of tree CSVs; supports multiple configs in one invocation (loads each tree once, runs hemo per config) |
| `tree_diagnostic.jl` | per-tree root R, path R statistics, diameter histogram |
| `dead_segments_diag.jl` | CSV (parent-id) vs FlowTree (BFS) reachability + vertex merge histogram — used to find the original load_tree bug; useful for any new CSV producer |
| `lcx_dead_flow_diag.jl` | per-segment flow distribution + top bottleneck list; useful when a tree gives unexpectedly low Q |
| `xcat_label_dims.jl` | groups segments by `label` to spot anatomical kinks (e.g. an XCAT chain that narrows mid-way) |
| `lad_stenosis_sweep.jl` + `plot_lad_stenosis.py` | proximal-LAD stenosis vs root flow at fixed dilation state (Gould curve); accepts a comma-separated stenosis list via `ARGS[4]` |
| `lad_stenosis_autoreg_sweep.jl` + `plot_autoreg_dilation.py` | closed-loop autoregulation: solves analytically for the arteriolar dilation factor that maintains target Q as stenosis grows, up to the 1.6× reserve cap |
| `extract_peak_iodine.jl` | run hemo + simulate_contrast for each tree, find global peak time, write per-segment peak iodine + per-segment `arrival_time.f32` (consumed by `simple_dynamic_viewer.jl` and the perfusion-pipeline voxelizer) |
| `simple_dynamic_viewer.jl` | standalone dynamic-contrast HTML viewer at the ≥500 μm "main artery" level. Reads `arrival_time.f32` if present (physical hemo arrival); falls back to a Murray-velocity v(s) = v_root · D_s/D_root model. Does **not** use FlowContrastSim's runtime so it scales to the 357 M-segment phantom-grown trees that would OOM `build_contrast_viewer` |
| `build_dynamic_contrast_viewer.jl` | full-runtime viewer (calls `build_contrast_viewer`); reduced budgets via `max_segments_per_branch = 800` + `time_stride = 10`. OOMs on the very largest 110 M-segment trees; use `simple_dynamic_viewer.jl` instead in that case |
| `run_8um_uniphase.jl` | full 8 μm phantom-grown tree run (~100 M segs / tree) using `coronary_hyperemic_uniphase.toml` + Taylor-Aris PDE + peak-only V2 snapshot. Writes `output_8um_uniphase/summary.txt` with per-tree root flow, arrival times, Taylor σ statistics, and V2 peak iodine concentration. Sparse mode at 50 μm keeps memory bounded |

All diagnostic scripts take the tree dir as `ARGS[1]` and never call
`run_flow_simulation` — they're fast (~5 min) compared to the full pipeline.

---

## Reproducing the canonical vmale50 coronary run

```bash
# step 1: grow trees with VascularTreeSim (~6 h)
cd ../VascularTreeSim.jl
julia --project=. --threads=auto examples/run_coronary_growth.jl configs/coronary.toml

# step 2: scale to at-rest (~7 min)
for tree in lad lcx rca; do
  julia --project=. scripts/scale_to_rest.jl \
        output/${tree}_segments.csv \
        output_at_rest/${tree}_segments.csv \
        0.55
done

# step 3: flow report for both states (~30 min)
cd ../FlowContrastSim.jl
julia --project=. scripts/natural_flow_summary.jl \
      ../VascularTreeSim.jl/output_at_rest configs/coronary_baseline.toml
julia --project=. scripts/natural_flow_summary.jl \
      ../VascularTreeSim.jl/output configs/coronary_hyperemic.toml
```

Expected per-tree numbers (Wong-Molloi 2008 alignment: scale_to_rest tone
= 0.375 band [8, 400] μm = 1.6× reserve; baseline cap_R = 0.15 + hyperemic
cap_R = 0.12):

| tree | baseline (mL/min) | hyperemic (mL/min) | CFR |
|---|---|---|---|
| LAD | 58.4 | 184.0 | 3.15× |
| LCX | 64.9 | 199.9 | 3.08× |
| RCA | 71.5 | 220.5 | 3.08× |
| **total** | **195** | **604** | — |

Per-tree baseline runs slightly above the 30-60 mL/min clinical band
because Wong-Molloi's empirical 1.6× was calibrated to Pantely 1984 /
Fearon 2004 swine data, not our exact tree resistance; raising the
factor (in scale_to_rest) trades off the CT-visible diameter realism.
Hyperemic 184-220 sits in the 120-240 clinical band, CFR 3.0-3.2× in the
3-5× literature range. Total flow ~195 (rest) / 604 (max) matches
population means.

### LAD proximal stenosis sweep — autoregulation curve

`scripts/lad_stenosis_autoreg_sweep.jl` then sweeps proximal LAD stenosis
in 10 % increments (plus 82, 84, 86, 88, 92, 94, 96, 98 % to resolve the
reserve-exhaustion knee). At each stenosis level the in-band (8-400 μm)
arterioles dilate up to 1.6× to maintain baseline Q ≈ 58.4 mL/min, until
reserve runs out:

| stenosis | baseline | dilation used | exhausted | hyperemic |
|---|---|---|---|---|
| 0-70 %  | 58.4 | 0-6 % | no | 184 → 135 |
| 80 %    | 58.7 | 40 % | no, close to limit | 70 |
| **82 %** | **52.8** | **60 %** | **YES — first** | 55 |
| 90 %    | 14.9 | 60 % (saturated) | yes | 16 |
| 98 %    | 0.14 | 60 % (saturated) | yes | 0.15 |

The plot at `scripts/lad_stenosis.png` shows the classic autoregulation +
Gould curve: baseline stays flat through 80 % stenosis, then drops sharply
at 82 % when the 1.6× arteriolar reserve is exhausted. Hyperemic (no
autoregulation, max-dilated geometry throughout) tracks the Gould curve
falling from ~50 % stenosis. The two curves merge at ≥ 82 % stenosis
because both run f = 1.0 once autoregulation can no longer help.

---

## References

- Murray CD. "The physiological principle of minimum work…" PNAS 1926
- Pries AR et al. "Blood viscosity in tube flow: dependence on diameter
  and hematocrit." Am J Physiol 1992;263:H1770-8
- Pries AR, Secomb TW. "Microvascular blood viscosity in vivo and the
  endothelial surface layer." Am J Physiol 2005;289:H2657-64
- Dodge JT et al. "Lumen diameter of normal human coronary arteries:
  influence of age, sex, anatomic variation, and left ventricular
  hypertrophy." Circulation 1992;86:232-46
- Gould KL. "Pressure-flow characteristics of coronary stenoses in
  unsedated dogs at rest and during coronary vasodilation." Circ Res 1978
- Molloi S, Wong JT. "Regional blood flow analysis and its relationship
  with arterial branch lengths and lumen volume in the coronary arterial
  tree." Phys Med Biol 2007;52:1495 — stem flow vs crown length/volume
  scaling, Mittal-Kassab tree reconstruction down to 8 μm.
- Wong JT, Molloi S. "Determination of fractional flow reserve (FFR)
  based on scaling laws: a simulation study." Phys Med Biol
  2008;53:3995-4011 — empirical 1.6× arteriolar (≤ 400 μm) dilation
  factor used here.
- Cornelissen AJM et al. "Myogenic reactivity and resistance distribution
  in the coronary arterial tree…" Am J Physiol 2000;278:H1490-9 —
  passive pressure-diameter curve (eq. 12 in Wong-Molloi 2008); referenced
  in the at-rest scaling but not yet implemented as iterative coupling
  with the flow solver.

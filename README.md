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
│   cap_R = 0.15          │                  │   cap_R = 0.15          │
│   (state diff comes     │                  │   (same literature      │
│    from tree geometry)  │                  │    value — no calib)    │
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
   baseline 47/52/57 mL/min               hyperemic 166/179/196 mL/min
   (LAD/LCX/RCA, CFR 3.4-3.5×)
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

**Strict dilation, not lumped tone.** Two coronary states use two trees:

- **Max-dilated**: the as-grown tree (Murray-optimal diameters). Used for the
  hyperemic config.
- **At-rest**: max-dilated tree with arteriolar diameters scaled down by the
  bell-curve tone function in
  [`VascularTreeSim.jl/scripts/scale_to_rest.jl`](../VascularTreeSim.jl/scripts/scale_to_rest.jl).
  Peak constriction at 100 μm, conduit (> 1 mm) and capillary (< 10 μm)
  almost unchanged.

Both configs use the same `cap_R = 0.15`. CFR (coronary flow reserve) is then
a property of the **anatomical scaling**, not a free parameter — matching how
real myocardium controls perfusion through arteriolar smooth muscle. This
also gives the right vessel diameters to a virtual CT scanner: an at-rest CT
sees constricted arterioles, a hyperemic CT sees max-dilated arterioles, both
from the same underlying anatomy.

**Sparse contrast.** The default `contrast_min_diameter_um = 50` skips
capillary-level segments from the `(n_segs × n_timesteps)` concentration
matrix. Without this a 110 M-segment LCX would need ~50 GB just for the
contrast field. Hemodynamics (R + flow + transit time) is still computed for
every segment.

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

---

## Config TOML schema

```toml
root_pressure_mmhg    = 100.0
terminal_pressure_mmhg = 15.0
discharge_hematocrit  = 0.45

# physics-grounded downstream R (Pries-Secomb in-vivo for myocardium):
capillary_bed_R_per_100g_mmHgmin_ml = 0.15

# contrast bolus shape (gamma-variate)
contrast_amplitude = 5.0   # mg I / mL peak
contrast_t0   = 0.5        # s onset
contrast_tmax = 4.0        # s peak
contrast_alpha = 3.0       # shape
dt    = 0.1                # s — simulation timestep
t_end = 20.0               # s
max_arrival_s = 15.0       # truncate plug-flow arrival lookup
contrast_min_diameter_um = 50.0   # skip contrast on segs below this diameter

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
```

Both `configs/coronary_baseline.toml` and `configs/coronary_hyperemic.toml`
ship with the same `cap_R = 0.15` and `territory_masses_g`; only the comment
about which tree directory to pair them with differs. To target a brain
tree, swap the masses, swap the tree dir, leave the physics knobs alone.

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

`contrast_results[name]` is a `ContrastResult` with a sparse `(n_sim_segs ×
n_timesteps)` concentration matrix and a `simulated_segment_ids` index.
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
| `lad_stenosis_sweep.jl` + `plot_lad_stenosis.py` | proximal-LAD stenosis vs root flow (Gould curve) |

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

Expected (same numbers as the README's other repo):

| tree | baseline (mL/min) | hyperemic (mL/min) | CFR |
|---|---|---|---|
| LAD | 47 | 166 | 3.5× |
| LCX | 52 | 179 | 3.4× |
| RCA | 57 | 196 | 3.4× |
| total | **157** | **541** | — |

All three trees fall in the clinical bands (baseline 30–60 mL/min,
hyperemic 120–240 mL/min, CFR 3–5×). Total flow matches population means
(~150 mL/min at rest, 400–700 mL/min peak hyperemia).

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

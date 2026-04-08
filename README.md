# FlowContrastSim.jl

A Julia package for hemodynamic flow simulation and dynamic contrast transport
in vascular trees.

Given a set of vascular segment CSVs produced by the companion
**VascularTreeSim.jl** package (or any CSV with the same schema), this package
computes Poiseuille-based blood flow with diameter-dependent non-Newtonian
viscosity (Pries 1992 / Pries & Secomb 2005), calibrates terminal capillary-bed
resistance to match a target flow rate, and simulates a bolus of iodinated
contrast agent propagating through the network. Results are exported to a
self-contained interactive 3-D HTML viewer.

---

## Features

- **Poiseuille hemodynamics** -- Hagen-Poiseuille resistance for every segment, with parallel-resistance reduction at bifurcations and top-down pressure/flow computation.
- **Pries non-Newtonian viscosity** -- Diameter- and hematocrit-dependent apparent viscosity using the Pries 1992 in-vitro model plus the Pries & Secomb 2005 endothelial surface layer (ESL) correction, capturing the Fahraeus-Lindqvist effect in microvessels.
- **Terminal capillary-bed calibration** -- A lumped resistance is added to every terminal segment and tuned via binary search so that total tree inflow matches a user-specified target flow (e.g. 242 mL/min for the LAD under hyperemic conditions).
- **Gamma-variate contrast injection** -- Configurable bolus shape (amplitude, onset, peak time, shape exponent).
- **Plug-flow contrast transport** -- Arrival-time computation with log-compressed time mapping and Taylor-like dispersion broadening for deeper segments.
- **Interactive 3-D HTML viewer** -- Plotly.js-based visualization with time slider, play/pause animation, per-segment hover info (branch, diameter, length), and a heat-map color scale showing local contrast concentration.

---

## Installation

Requires Julia >= 1.9. Not yet registered in the Julia General registry.

```julia
using Pkg

# From a local clone:
Pkg.develop(path="/path/to/FlowContrastSim.jl")

# Or from a Git URL:
Pkg.add(url="https://github.com/your-org/FlowContrastSim.jl.git")
```

## Testing

```julia
using Pkg
Pkg.test("FlowContrastSim")
```

Runs unit tests for the Pries viscosity model, gamma-variate input, synthetic
tree hemodynamics (Y-bifurcation flow conservation), contrast transport, config
loading, and viewer generation.

### Dependencies

All dependencies are standard Julia packages declared in `Project.toml`:

| Package | Purpose |
|---|---|
| `DelimitedFiles` | CSV reading |
| `LinearAlgebra` | Vector norms for segment lengths |
| `StaticArrays` | `SVector{3,Float64}` for 3-D vertex coordinates |
| `Statistics` | Summary statistics |
| `TOML` | Configuration file parsing |

No external binaries or compiled libraries are required.

---

## Quick Start

```julia
using FlowContrastSim

# 1. Load configuration
config = load_flow_config("configs/coronary_hyperemic.toml")

# 2. Run the full pipeline (discover CSVs, compute flow, simulate contrast, build viewer)
results = run_flow_simulation("path/to/tree_csvs/", config; output_dir="output")

# results.trees             -- Dict{String, FlowTree}
# results.hemo_results      -- Dict{String, HemodynamicsResult}
# results.contrast_results  -- Dict{String, ContrastResult}
# results.viewer_path       -- path to the generated HTML file
```

Or drive each step manually:

```julia
using FlowContrastSim

# Load a single tree
tree = load_tree("LAD", "lad_grown_segments.csv")

# Compute hemodynamics with a target flow of 242 mL/min
hemo = compute_hemodynamics(tree;
    root_pressure   = 100.0 * 133.322,   # 100 mmHg in Pa
    terminal_pressure = 15.0 * 133.322,   # 15 mmHg in Pa
    hematocrit      = 0.45,
    target_flow_ml_min = 242.0)

# Simulate contrast bolus
cr = simulate_contrast(tree, hemo;
    dt=0.1, t_end=20.0,
    amplitude=5.0, t0=0.5, tmax=4.0, alpha=3.0)

# Access per-segment data
hemo.segment_flow          # m^3/s per segment
hemo.pressure_proximal     # Pa
hemo.transit_time_s        # seconds
cr.concentration           # [n_segments x n_timesteps] matrix, mg/mL
```

---

## Input Format

Tree CSVs are produced by the **VascularTreeSim.jl** package. Each row describes one vessel segment. Required columns:

| Column | Units | Description |
|---|---|---|
| `segment_id` | integer | Unique segment identifier |
| `parent_segment_id` | integer | ID of the parent segment (`-1` or empty for root segments) |
| `x1_cm`, `y1_cm`, `z1_cm` | cm | Proximal (start) endpoint coordinates |
| `x2_cm`, `y2_cm`, `z2_cm` | cm | Distal (end) endpoint coordinates |
| `diameter_um` | micrometers | Segment diameter |
| `label` | string | Segment category (e.g. `"epicardial"`, `"grown"`) |

Additional columns present in the CSV but not consumed by this package: `branch`, `xmid_cm`, `ymid_cm`, `zmid_cm`, `length_mm`.

When `parent_segment_id` is present, topology is reconstructed deterministically via the parent-child links. When it is absent (older CSVs), the loader falls back to nearest-vertex matching for backward compatibility.

---

## FlowConfig TOML Schema

All simulation parameters are set in a single TOML file. Below is the complete schema with default values:

```toml
# ── Hemodynamics ──
root_pressure_mmhg     = 100.0   # Aortic root pressure (mmHg)
terminal_pressure_mmhg = 15.0    # Distal capillary bed pressure (mmHg)
discharge_hematocrit   = 0.45    # Systemic discharge hematocrit (fraction, 0-1)

# ── Time grid ──
dt    = 0.1    # Simulation time step (seconds)
t_end = 20.0   # Total simulation duration (seconds)

# ── Contrast bolus (gamma-variate) ──
contrast_amplitude = 5.0    # Peak concentration at injection site (mg/mL)
contrast_t0        = 0.5    # Bolus onset time (seconds)
contrast_tmax      = 4.0    # Time of peak input concentration (seconds)
contrast_alpha     = 3.0    # Shape exponent (higher = sharper bolus)

# ── Contrast transport ──
max_arrival_s = 15.0   # Maximum compressed arrival time window (seconds).
                       # Controls the time-compression mapping that maps
                       # raw arrival times into the visualization window.
                       # Set to 0 for automatic (t_end - tmax).

# ── Target flow rates (per branch) ──
# Keys must match branch names derived from CSV filenames.
# Branch names are uppercased, with suffixes like "_grown_segments"
# stripped automatically (e.g. "lad_grown_segments.csv" -> "LAD").
# Set to 0 or omit a branch to skip calibration for that branch.
[target_flows_ml_min]
LAD = 242.0    # mL/min, hyperemic LAD flow
LCX = 116.0    # mL/min, hyperemic LCx flow
RCA = 214.0    # mL/min, hyperemic RCA flow
```

### Field Reference

| Field | Type | Default | Description |
|---|---|---|---|
| `root_pressure_mmhg` | Float | 100.0 | Inlet pressure at the root of the tree |
| `terminal_pressure_mmhg` | Float | 15.0 | Outlet pressure at every terminal (capillary bed) |
| `discharge_hematocrit` | Float | 0.45 | Hematocrit used for Pries viscosity calculation |
| `dt` | Float | 0.1 | Contrast simulation time step in seconds |
| `t_end` | Float | 20.0 | End time of contrast simulation in seconds |
| `contrast_amplitude` | Float | 5.0 | Peak bolus concentration (mg/mL) |
| `contrast_t0` | Float | 0.5 | Bolus arrival time at tree root (seconds) |
| `contrast_tmax` | Float | 4.0 | Time of peak input concentration (seconds) |
| `contrast_alpha` | Float | 3.0 | Gamma-variate shape parameter |
| `max_arrival_s` | Float | 15.0 | Compressed arrival-time window (seconds); 0 = auto |
| `target_flows_ml_min` | Table | `{}` | Per-branch target total flow in mL/min |

---

## Hemodynamics Model

### Overview

The hemodynamic solver treats the vascular tree as an electrical circuit: each vessel segment is a resistor, bifurcations are parallel junctions, and flow is driven by the pressure difference between the aortic root and the terminal capillary bed.

### Poiseuille Resistance

Each segment's resistance is:

```
R = 8 * mu * L / (pi * r^4)
```

where `L` is the segment length, `r` is the radius, and `mu` is the **apparent viscosity** (not a fixed constant -- see below).

### Diameter-Dependent Viscosity (Pries Model)

Blood is not a Newtonian fluid. In tubes smaller than about 300 micrometers, the apparent viscosity drops significantly due to the formation of a cell-free plasma layer near the vessel wall (the **Fahraeus-Lindqvist effect**). Below roughly 7 micrometers, viscosity rises sharply again because red blood cells must deform to squeeze through. This package implements the full Pries parameterization rather than using a constant bulk viscosity.

The computation proceeds in three steps:

1. **Reference viscosity at Hd = 0.45** (Pries 1992, Equations 1-3):

   ```
   eta_0.45(D) = 220 * exp(-1.3*D) + 3.2 - 2.44 * exp(-0.06 * D^0.645)
   ```

   where `D` is the vessel diameter in micrometers. This captures the U-shaped viscosity curve: high for very small tubes, minimum around 7 micrometers, rising back toward the bulk value for large tubes.

2. **Hematocrit correction** via the `C(D)` parameter:

   ```
   C(D) = (0.8 + exp(-0.075*D)) * (-1 + 1/(1 + 1e-11 * D^12)) + 1/(1 + 1e-11 * D^12)

   eta_rel(D, Hd) = 1 + (eta_0.45 - 1) * ((1 - Hd)^C - 1) / ((1 - 0.45)^C - 1)
   ```

3. **Endothelial surface layer (ESL) correction** (Pries & Secomb 2005): For vessels with diameters between 10 and 150 micrometers, a 1.1 micrometer glycocalyx layer lines the endothelium, reducing the effective lumen:

   ```
   D_eff = D - 2 * 1.1 um
   eta_in_vivo = eta_in_vitro * (D / D_eff)^4
   ```

   The fourth-power scaling comes from the Poiseuille relation: halving the effective radius increases resistance 16-fold.

The apparent viscosity in Pa*s is then:

```
mu = eta_plasma * eta_rel
```

where `eta_plasma = 0.0012 Pa*s` (1.2 cP).

For vessels below 2.7 micrometers, the model returns an effectively infinite viscosity (flow is blocked).

### Subtree Resistance and Flow Distribution

Resistances are combined bottom-up:

- **Terminal segments**: `R_subtree = R_segment + R_bed` (where `R_bed` is the calibrated capillary-bed resistance).
- **Bifurcations**: Children combine in parallel: `1/R_children = sum(1/R_child_i)`, then `R_subtree = R_segment + R_children`.

Flow is distributed top-down: the root segment receives `Q = dP / R_subtree_root`, and at each bifurcation flow splits inversely proportional to child subtree resistances.

### Terminal Bed Calibration

When a `target_flow_ml_min` is specified for a branch, the solver finds a uniform terminal-bed resistance `R_bed` (appended to every leaf segment) such that total tree inflow equals the target. This is solved by binary search over 60 iterations, which converges to machine precision. If the tree's intrinsic resistance already limits flow below the target, `R_bed` is set to zero.

---

## Contrast Transport Model

### Plug-Flow Assumption

Contrast agent is modeled as a passive tracer that propagates with the local blood velocity. Each segment's **transit time** is:

```
tau = V / Q = (pi * r^2 * L) / Q
```

The **arrival time** at a segment's midpoint is the cumulative sum of half-transit-times along the path from root:

```
arrival[root] = tau_root / 2
arrival[child] = arrival[parent] + tau_parent/2 + tau_child/2
```

### Time Compression

Raw arrival times in a large tree can span many orders of magnitude (microseconds in epicardial vessels to minutes in terminal arterioles). For visualization, arrival times are log-compressed into a configurable window:

```
arrival_compressed = max_arrival_s * log(1 + arrival_raw / scale) / log(1 + p99 / scale)
```

where `scale` is the 25th-percentile raw arrival time and `p99` is the 99th percentile.

### Input Bolus

The contrast input at the tree root follows a **gamma-variate** curve:

```
C_input(t) = A * ((t - t0) / (tmax - t0))^alpha * exp(alpha * (1 - (t - t0) / (tmax - t0)))
```

for `t > t0`, and zero otherwise. Parameters:
- `A` (amplitude): peak concentration in mg/mL
- `t0`: bolus onset time
- `tmax`: time of peak concentration
- `alpha`: shape exponent (higher values produce a sharper, more peaked bolus)

### Dispersion

As contrast propagates deeper into the tree, the bolus broadens due to Taylor dispersion and mixing at bifurcations. This is modeled by a dispersion factor:

```
disp_factor = sqrt(1 + arrival / t_disp)
```

where `t_disp = 3.0 s` is the dispersion time constant. The factor stretches the time axis of the input curve and reduces its amplitude by `1/disp_factor`, conserving the total amount of contrast.

---

## Viewer

`build_contrast_viewer` generates a self-contained HTML file that uses Plotly.js (loaded from CDN) for 3-D rendering. No server is needed; open the file in any modern browser.

### Features

- **Time slider** -- Drag to any time point in the simulation. The readout shows the current time in seconds.
- **Play / Pause** -- Animate through all time frames at approximately 12 fps (80 ms per frame).
- **3-D rotation, zoom, pan** -- Standard Plotly.js orbit controls (click-drag to rotate, scroll to zoom, right-drag to pan).
- **Segment hover info** -- Hovering over a marker shows: branch name, segment ID, diameter in micrometers, and length in millimeters.
- **Concentration color scale** -- Segments are colored from dark (no contrast) through blue and orange to red (peak concentration). A color bar shows the mapping in mg/mL.
- **Skeleton overlay** -- A faint wireframe of the full tree structure is drawn beneath the concentration markers so that unfilled segments remain visible.
- **Multi-branch support** -- Each branch (LAD, LCx, RCA, etc.) is rendered as a separate Plotly trace with an auto-assigned color from a 10-color palette. Branches can be toggled on and off via the legend.
- **Segment budget** -- To keep the viewer responsive, each branch is limited to `max_segments_per_branch` marker points (default 6000). Non-grown (epicardial) segments are always included; grown segments are prioritized by diameter. A subsampled skeleton is added for visual completeness.

---

## API Reference

### Types

**`FlowTree`** -- Read-only representation of a vascular tree loaded from CSV.

| Field | Type | Description |
|---|---|---|
| `name` | `String` | Branch name |
| `vertices` | `Vector{SVector{3,Float64}}` | 3-D vertex coordinates in cm |
| `parent_vertex` | `Vector{Int}` | Parent vertex index for each vertex |
| `children` | `Vector{Vector{Int}}` | Child vertex indices for each vertex |
| `incoming_segment` | `Vector{Int}` | Segment index whose distal end is this vertex |
| `segment_start` | `Vector{Int}` | Proximal vertex index per segment |
| `segment_end` | `Vector{Int}` | Distal vertex index per segment |
| `segment_diameter_cm` | `Vector{Float64}` | Diameter in cm |
| `segment_label` | `Vector{String}` | Label string (e.g. `"epicardial"`, `"grown"`) |
| `root_vertex` | `Int` | Index of the root vertex |

**`FlowConfig`** -- All simulation parameters (see TOML schema above).

**`HemodynamicsResult`** -- Output of `compute_hemodynamics`.

| Field | Type | Units | Description |
|---|---|---|---|
| `segment_resistance` | `Vector{Float64}` | Pa*s/m^3 | Poiseuille resistance per segment |
| `segment_flow` | `Vector{Float64}` | m^3/s | Volume flow rate per segment |
| `pressure_proximal` | `Vector{Float64}` | Pa | Pressure at proximal end |
| `pressure_distal` | `Vector{Float64}` | Pa | Pressure at distal end |
| `segment_volume_m3` | `Vector{Float64}` | m^3 | Cylindrical volume of each segment |
| `transit_time_s` | `Vector{Float64}` | s | Volume / flow rate |

**`ContrastResult`** -- Output of `simulate_contrast`.

| Field | Type | Description |
|---|---|---|
| `times` | `Vector{Float64}` | Time grid in seconds |
| `concentration` | `Matrix{Float64}` | `[n_segments x n_timesteps]` concentration in mg/mL |
| `outlet_concentration` | `Matrix{Float64}` | `[n_segments x n_timesteps]` outlet concentration in mg/mL |

### Functions

```julia
load_flow_config(path::String) -> FlowConfig
```
Parse a TOML configuration file into a `FlowConfig`.

```julia
load_tree(name::String, csv_path::String) -> FlowTree
```
Load a single vascular tree from a CSV file. `name` is an arbitrary label (e.g. `"LAD"`).

```julia
load_trees(dict::Dict{String,String}) -> Dict{String, FlowTree}
```
Load multiple trees from a dictionary of `name => csv_path`.

```julia
compute_hemodynamics(tree::FlowTree;
    root_pressure::Float64       = 13332.0,   # 100 mmHg in Pa
    terminal_pressure::Float64   = 1999.8,     # 15 mmHg in Pa
    hematocrit::Float64          = 0.45,
    target_flow_ml_min::Float64  = 0.0         # 0 = no calibration
) -> HemodynamicsResult
```
Compute steady-state Poiseuille flow with Pries viscosity. If `target_flow_ml_min > 0`, a terminal-bed resistance is calibrated so that total root inflow matches the target.

```julia
simulate_contrast(tree::FlowTree, hemo::HemodynamicsResult;
    dt::Float64            = 0.05,
    t_end::Float64         = 30.0,
    root_input             = nothing,    # custom input curve, or nothing for gamma-variate
    amplitude::Float64     = 5.0,
    t0::Float64            = 0.5,
    tmax::Float64          = 4.0,
    alpha::Float64         = 3.0,
    max_arrival_s::Float64 = 0.0         # 0 = auto
) -> ContrastResult
```
Simulate contrast propagation using plug-flow transport with dispersion. Pass a custom `root_input` vector (same length as `times`) to override the default gamma-variate bolus.

```julia
gamma_variate_input(times::Vector{Float64};
    amplitude=5.0, t0=0.5, tmax=4.0, alpha=3.0
) -> Vector{Float64}
```
Generate a gamma-variate concentration-time curve. Useful for constructing a custom input or for plotting the injection profile.

```julia
build_contrast_viewer(path, trees, hemo_results, contrast_results;
    time_stride=1, title="Dynamic Contrast Transport",
    max_segments_per_branch=8000, branch_colors=nothing
) -> String
```
Write a self-contained HTML file with an interactive 3-D contrast viewer. Returns the output path. `time_stride` controls temporal subsampling (e.g. `3` keeps every third frame to reduce file size). `branch_colors` is an optional `Dict{String,String}` mapping branch names to CSS color strings.

```julia
run_flow_simulation(tree_csv_dir::String, config::FlowConfig;
    output_dir="output"
) -> NamedTuple
```
End-to-end pipeline: discover CSVs in `tree_csv_dir`, compute hemodynamics, simulate contrast, build viewer, and return all results. CSV files are matched by `*.csv`; branch names are derived from filenames by stripping suffixes like `_grown_segments` and uppercasing.

```julia
pries_viscosity_relative(diameter_um::Float64; hematocrit=0.45) -> Float64
```
Compute the relative apparent viscosity (eta_apparent / eta_plasma) for a given vessel diameter and hematocrit using the Pries 1992 model with the Pries & Secomb 2005 ESL correction.

```julia
apparent_viscosity(diameter_um::Float64; hematocrit=0.45) -> Float64
```
Compute the absolute apparent viscosity in Pa*s (= eta_plasma * eta_relative).

---

## Physics: Why Constant Viscosity Is Wrong for Microvessels

### The Fahraeus-Lindqvist Effect

When blood flows through tubes smaller than about 300 micrometers in diameter, its apparent viscosity is significantly lower than the bulk value measured in a standard viscometer (approximately 3.5 cP at 45% hematocrit). This phenomenon, first reported by Fahraeus and Lindqvist in 1931, arises because red blood cells migrate toward the tube axis, leaving a cell-free (or cell-depleted) plasma layer along the wall. Since the wall shear rate is highest near the boundary, and the near-wall fluid is predominantly low-viscosity plasma, the effective (apparent) viscosity of the tube as a whole drops.

The effect is diameter-dependent:

- **D > 300 um**: The cell-free layer is thin relative to the tube diameter. Apparent viscosity approaches the bulk value.
- **D ~ 10-50 um**: The cell-free layer occupies a significant fraction of the cross-section. Apparent viscosity can fall to less than half the bulk value.
- **D ~ 5-7 um**: Viscosity reaches a minimum. Below this, red blood cells (diameter ~6-8 um but highly deformable) must squeeze single-file through the lumen, and viscosity rises steeply.
- **D < 2.7 um**: Even maximally deformed RBCs cannot transit. Flow is effectively impossible.

### The Pries 1992 Model

Pries, Neuhaus, and Gaehtgens (1992) compiled measurements of blood viscosity in glass tubes across a wide range of diameters (3.3-1978 um) and hematocrits (0-0.90). They fit a parametric model with three components:

1. A reference curve `eta_0.45(D)` for blood at 45% hematocrit.
2. A diameter-dependent exponent `C(D)` that governs how viscosity scales with hematocrit.
3. The full viscosity `eta(D, Hd)` interpolating between plasma (Hd=0) and the reference curve.

This in-vitro model accurately captures the Fahraeus-Lindqvist minimum and the steep rise at small diameters.

### The Endothelial Surface Layer (Pries & Secomb 2005)

Glass tubes lack the glycocalyx that lines living endothelium. In vivo, this endothelial surface layer (ESL) extends approximately 1.1 micrometers into the lumen, reducing the effective diameter available for flow. For a 20-micrometer vessel, this means the effective diameter is 17.8 micrometers -- and because Poiseuille resistance scales with the fourth power of radius, this raises resistance by `(20/17.8)^4 = 1.59` (a 59% increase).

The ESL correction is applied for diameters between 10 and 150 micrometers. Below 10 micrometers, the glycocalyx is compressed by transiting red blood cells and has minimal additional effect. Above 150 micrometers, the 1.1-micrometer layer is negligible relative to the lumen.

### Why This Matters

Using a constant viscosity of 3.5 cP for all vessel sizes produces systematic errors:

- **Overestimates resistance in 10-100 um vessels** by up to a factor of 2, because it ignores the Fahraeus-Lindqvist reduction.
- **Underestimates resistance in 3-7 um vessels**, where single-file RBC transit dramatically increases effective viscosity.
- **Ignores the ESL**, which adds 20-60% additional resistance in vessels of 10-50 micrometers.

In a tree with thousands of microvessels, these errors compound and can shift the total flow by 30-50%, alter the flow distribution at bifurcations, and produce incorrect transit times for contrast simulation.

---

## References

1. **Pries AR, Neuhaus D, Gaehtgens P.** Blood viscosity in tube flow: dependence on diameter and hematocrit. *American Journal of Physiology -- Heart and Circulatory Physiology*. 1992;263(6):H1770-H1778. -- The foundational parameterization of in-vitro apparent blood viscosity as a function of tube diameter and hematocrit.

2. **Pries AR, Secomb TW.** Microvascular blood viscosity in vivo and the endothelial surface layer. *American Journal of Physiology -- Heart and Circulatory Physiology*. 2005;289(6):H2657-H2664. -- Extension of the 1992 model to account for the endothelial glycocalyx layer that reduces effective lumen diameter in vivo.

3. **Molloi S, et al.** Estimation of coronary artery hyperemic blood flow based on arterial lumen volume using angiographic images. *International Journal of Cardiovascular Imaging*. 2007. -- Basis for the target hyperemic flow rates used in calibration.

4. **Wong JT, et al.** Quantification of fractional flow reserve based on angiographic image data. *International Journal of Cardiovascular Imaging*. 2008. -- Related work on flow quantification from angiographic data.

5. **Fahraeus R, Lindqvist T.** The viscosity of the blood in narrow capillary tubes. *American Journal of Physiology*. 1931;96(3):562-568. -- Original observation of the diameter-dependent viscosity reduction in narrow tubes.

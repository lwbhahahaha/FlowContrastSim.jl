#!/usr/bin/env python3
"""Plot LAD stenosis sweep: stenosis % vs root flow (mL/min).

Reads scripts/lad_stenosis.csv produced by scripts/lad_stenosis_sweep.jl.
Writes scripts/lad_stenosis.png.
"""

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
CSV_IN = HERE / "lad_stenosis.csv"
PNG_OUT = HERE / "lad_stenosis.png"

rows = []
with open(CSV_IN, newline="") as fh:
    reader = csv.DictReader(fh)
    for row in reader:
        rows.append(
            (
                float(row["stenosis_pct"]),
                float(row["proximal_diameter_um"]),
                float(row["root_flow_mlmin"]),
            )
        )

stenosis = [r[0] for r in rows]
flow = [r[2] for r in rows]
baseline = flow[0]
rel = [100.0 * q / baseline for q in flow]

fig, ax1 = plt.subplots(figsize=(9, 5.5), dpi=120)

# Absolute flow curve
color_abs = "#1f77b4"
ax1.plot(stenosis, flow, "-o", color=color_abs, linewidth=2.4, markersize=7, label="Root flow (mL/min)")
ax1.set_xlabel("Proximal LAD stenosis (%)", fontsize=12)
ax1.set_ylabel("Root flow (mL/min)", fontsize=12, color=color_abs)
ax1.tick_params(axis="y", labelcolor=color_abs)
ax1.set_xlim(-3, 93)
ax1.set_ylim(0, max(flow) * 1.08)
ax1.grid(True, alpha=0.3)
ax1.set_xticks(range(0, 100, 10))

# Relative flow on twin axis
ax2 = ax1.twinx()
color_rel = "#d62728"
ax2.plot(stenosis, rel, "--s", color=color_rel, linewidth=1.2, markersize=4.5, alpha=0.55, label="Relative (% baseline)")
ax2.set_ylabel("Relative flow (% of baseline)", fontsize=12, color=color_rel)
ax2.tick_params(axis="y", labelcolor=color_rel)
ax2.set_ylim(0, 110)

# Annotate baseline + clinically significant threshold
ax1.axhline(baseline, color=color_abs, linestyle=":", linewidth=0.9, alpha=0.5)
ax1.text(2, baseline + 6, f"baseline {baseline:.1f} mL/min", color=color_abs, fontsize=9)

ax1.axvspan(70, 90, color="#ffe5e5", alpha=0.4, zorder=0)
ax1.axvline(70, color="#b03030", linestyle="--", linewidth=1.0, alpha=0.6)
ax1.text(71, max(flow) * 0.92, "clinically significant\n(≥70% stenosis)", color="#8a2020", fontsize=9)

# Annotate a few knee points
for x, y in zip(stenosis, flow):
    if x in (0, 50, 70, 90):
        ax1.annotate(
            f"{y:.1f}",
            xy=(x, y),
            xytext=(6, 6),
            textcoords="offset points",
            fontsize=8,
            color=color_abs,
        )

plt.title("LAD stenosis sweep — natural forward hemodynamics\n(proximal dias_lad1 segments, 23 of them, uniform diameter reduction)", fontsize=11)
fig.tight_layout()
fig.savefig(PNG_OUT)
print(f"wrote {PNG_OUT}")

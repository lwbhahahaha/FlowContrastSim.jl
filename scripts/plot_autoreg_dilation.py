#!/usr/bin/env python3
"""Plot arteriolar dilation % vs proximal LAD stenosis %.

Reads scripts/lad_stenosis_baseline.csv (which has an autoreg_f column from
the autoregulation sweep) and writes a PNG showing how much of the 60 %
dilation reserve is recruited at each stenosis level.

  dilation_pct = (f / F_MIN − 1) × 100        with F_MIN = 0.625
       f = 0.625  → 0 %    (at-rest, full tone, reserve untouched)
       f = 1.000  → 60 %   (max dilation, reserve fully spent)

Usage:
    python plot_autoreg_dilation.py [baseline_csv] [output_png]
"""

import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

F_MIN = 0.625  # must match lad_stenosis_autoreg_sweep.jl

HERE = os.path.dirname(os.path.abspath(__file__))
CSV_IN  = sys.argv[1] if len(sys.argv) >= 2 else os.path.join(HERE, "lad_stenosis_baseline.csv")
PNG_OUT = sys.argv[2] if len(sys.argv) >= 3 else os.path.join(HERE, "lad_autoreg_dilation.png")

rows = []
with open(CSV_IN, newline="") as fh:
    reader = csv.DictReader(fh)
    if "autoreg_f" not in reader.fieldnames:
        sys.exit(f"{CSV_IN} has no autoreg_f column — was it produced by lad_stenosis_autoreg_sweep.jl?")
    for r in reader:
        rows.append((float(r["stenosis_pct"]),
                     float(r["autoreg_f"]),
                     int(r.get("reserve_exhausted", 0))))
rows.sort(key=lambda r: r[0])

stenosis = [r[0] for r in rows]
f_vals   = [r[1] for r in rows]
exhausted = [r[2] for r in rows]
dilation_pct = [(f / F_MIN - 1.0) * 100.0 for f in f_vals]

fig, ax = plt.subplots(figsize=(9, 5.5), dpi=120)

# Curve: dilation %
ax.plot(stenosis, dilation_pct, "-o", color="#4682b4", linewidth=2.4, markersize=7,
        label="Arteriolar dilation (% of reserve used)")

# Annotate every point with its dilation %
for x, y, ex in zip(stenosis, dilation_pct, exhausted):
    ax.annotate(f"{y:.1f}%", xy=(x, y), xytext=(5, 6),
                textcoords="offset points", fontsize=8, color="#264e6e")

# Mark the reserve-exhausted region
first_exhausted = next((s for s, e in zip(stenosis, exhausted) if e), None)
if first_exhausted is not None:
    ax.axvspan(first_exhausted, max(stenosis) + 3, color="#ffe5e5", alpha=0.4, zorder=0)
    ax.axvline(first_exhausted, color="#b03030", linestyle="--", linewidth=1.0, alpha=0.7)
    ax.text(first_exhausted + 0.5, 55, f"reserve exhausted\n(stenosis ≥ {first_exhausted:.0f}%)",
            color="#8a2020", fontsize=9, va="top")

# Reference lines: 0 % (at-rest) and 60 % (max reserve)
ax.axhline(0, color="#888", linestyle=":", linewidth=0.8, alpha=0.6)
ax.axhline(60, color="#888", linestyle=":", linewidth=0.8, alpha=0.6)
ax.text(-2, 1.5, "no dilation\n(at-rest)", fontsize=8, color="#666", va="bottom")
ax.text(-2, 58, "max reserve\n(60% dilation)", fontsize=8, color="#666", va="top")

ax.set_xlabel("Proximal LAD stenosis (%)", fontsize=12)
ax.set_ylabel("Arteriolar dilation (% of 60% reserve used)", fontsize=12)
ax.set_xlim(-3, max(stenosis) + 3)
ax.set_ylim(-3, 67)
ax.set_xticks(range(0, int(max(stenosis)) + 2, 10))
ax.grid(True, alpha=0.3)
ax.legend(fontsize=11, loc="upper left")

ax.set_title("Autoregulatory arteriolar dilation vs proximal LAD stenosis\n"
             "(8-400 μm band, 1.6× reserve cap)", fontsize=12)
fig.tight_layout()
os.makedirs(os.path.dirname(os.path.abspath(PNG_OUT)) or ".", exist_ok=True)
fig.savefig(PNG_OUT)
print(f"wrote {PNG_OUT}")

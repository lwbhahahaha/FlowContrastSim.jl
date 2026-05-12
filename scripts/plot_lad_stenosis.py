#!/usr/bin/env python3
"""Plot LAD proximal-stenosis sweep — flow vs % stenosis.

Reads one or two CSVs produced by `scripts/lad_stenosis_sweep.jl` and
emits a PNG. With two CSVs both curves are overlaid on shared axes
(baseline in blue, hyperemic in red) for direct CFR-vs-stenosis comparison.

Usage:
    # single state
    python plot_lad_stenosis.py scripts/lad_stenosis_coronary_baseline.csv

    # both states overlaid
    python plot_lad_stenosis.py \\
        scripts/lad_stenosis_coronary_baseline.csv \\
        scripts/lad_stenosis_coronary_hyperemic.csv \\
        scripts/lad_stenosis.png
"""

import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_csv(path):
    rows = []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            rows.append((float(row["stenosis_pct"]),
                         float(row["proximal_diameter_um"]),
                         float(row["root_flow_mlmin"])))
    rows.sort(key=lambda r: r[0])
    return rows


def label_from_path(path):
    base = os.path.splitext(os.path.basename(path))[0]
    # lad_stenosis_coronary_baseline -> "baseline"
    parts = base.split("_")
    if "baseline" in parts:
        return "baseline"
    if "hyperemic" in parts or "hyperemia" in parts:
        return "hyperemic"
    return base.replace("lad_stenosis_", "")


STATE_COLOR = {"baseline": "#1f77b4", "hyperemic": "#d62728"}


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    csv_paths = []
    out_png = None
    for a in sys.argv[1:]:
        if a.endswith(".csv"):
            csv_paths.append(a)
        elif a.endswith(".png"):
            out_png = a
    if not csv_paths:
        sys.exit("No CSV inputs.")
    if out_png is None:
        out_png = os.path.splitext(csv_paths[0])[0] + ".png"

    series = []
    for p in csv_paths:
        rows = load_csv(p)
        if not rows:
            print(f"warning: empty {p}")
            continue
        state = label_from_path(p)
        series.append((state, rows))

    fig, ax = plt.subplots(figsize=(9, 5.5), dpi=120)
    max_q = 0.0
    for state, rows in series:
        sten = [r[0] for r in rows]
        flow = [r[2] for r in rows]
        max_q = max(max_q, max(flow))
        color = STATE_COLOR.get(state, "#666666")
        ax.plot(sten, flow, "-o", color=color, linewidth=2.4, markersize=7, label=f"{state}")
        for x, y in zip(sten, flow):
            if x in (0, 30, 50, 70, 90):
                ax.annotate(f"{y:.1f}", xy=(x, y), xytext=(6, 6),
                            textcoords="offset points", fontsize=8, color=color)

    ax.set_xlabel("Proximal LAD stenosis (%)", fontsize=12)
    ax.set_ylabel("LAD root flow (mL/min)", fontsize=12)
    ax.set_xlim(-3, 93)
    ax.set_ylim(0, max_q * 1.10)
    ax.set_xticks(range(0, 100, 10))
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=11, loc="upper right")

    # clinical-significant threshold
    ax.axvspan(70, 93, color="#ffe5e5", alpha=0.4, zorder=0)
    ax.axvline(70, color="#b03030", linestyle="--", linewidth=1.0, alpha=0.6)
    ax.text(71, max_q * 0.94, "clinically significant\n(≥70% stenosis)", color="#8a2020", fontsize=9)

    if len(series) >= 2:
        title = "LAD root flow vs proximal stenosis — baseline vs hyperemic"
    else:
        title = f"LAD root flow vs proximal stenosis — {series[0][0]}"
    ax.set_title(title, fontsize=12)
    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(out_png)) or ".", exist_ok=True)
    fig.savefig(out_png)
    print(f"wrote {out_png}")


if __name__ == "__main__":
    main()

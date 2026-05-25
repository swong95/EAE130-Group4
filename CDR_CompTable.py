"""
CDR Poster Presentation — Aircraft Comparison Table
Exports as PNG (raster) or SVG (vector/scalable).

Fill in the VALUES dict below with your data, then run the script.
"""

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np

# ──────────────────────────────────────────────
#  DATA  ← fill these in
# ──────────────────────────────────────────────
AIRCRAFT = ["F/A-XX", "F/A-18 E/F", "Rafale-M"]

PARAMETERS = [
    "MTOW [lbs]",
    "Payload [lbs]",
    "Thrust : Weight",
    "Max Range [nmi]",
    "Cost [$USD 2026]",
]

# Row order matches PARAMETERS; column order matches AIRCRAFT
VALUES = [
    # F/A-XX        F/A-18 E/F      Rafale-M
    ["52,000",         "66,000",       "54,000"],   # MTOW
    ["16,600",         "17,750",       "21,000"],   # Payload
    ["0.83",         "0.93",         "1.00"],     # T:W
    ["1,900",         "1,000",        "1,800"],    # Max Range
    ["$73 M",         "$67.4 M",      "$95.5 M"],  # Cost
]

# ──────────────────────────────────────────────
#  STYLE
# ──────────────────────────────────────────────
COL_COLORS   = ["#0F1A2E", "#1B2A4A", "#2E5090", "#4472C4"]   # header bg per col (param + 3 a/c)
HEADER_FG    = "white"
ROW_COLORS   = ["#EEF2FA", "#FFFFFF"]   # alternating row fill
CELL_FG      = "#1A1A2E"
PARAM_FG     = "#1B2A4A"
FONT_FAMILY  = "DejaVu Sans"

OUTPUT_FILE  = r"d:\UC Davis\EAE 130AB\CDR_CompTable"   # extension added below
OUTPUT_FMT   = "png"    # "png"  or  "svg"
DPI          = 300      # only used for PNG

# ──────────────────────────────────────────────
#  BUILD TABLE
# ──────────────────────────────────────────────
n_rows = len(PARAMETERS)
n_cols = len(AIRCRAFT) + 1   # +1 for parameter label column

fig_w = 2.4 + 2.6 * len(AIRCRAFT)   # scale width with number of columns
fig_h = 0.70 * (n_rows + 1)         # tight: no title headroom
fig, ax = plt.subplots(figsize=(fig_w, fig_h))
ax.set_axis_off()

col_widths = [0.32] + [0.68 / len(AIRCRAFT)] * len(AIRCRAFT)   # relative fractions
col_x = np.cumsum([0.0] + col_widths[:-1])   # left edge of each column

ROW_H = 1.0 / (n_rows + 1)   # normalised height per row (header = row 0)

def draw_cell(ax, x, y, w, h, text, bg, fg, bold=False, fontsize=13):
    rect = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="square,pad=0",
        linewidth=0.8,
        edgecolor="#8A9CC2",
        facecolor=bg,
        transform=ax.transAxes,
        clip_on=False,
    )
    ax.add_patch(rect)
    ax.text(
        x + w / 2, y + h / 2, text,
        ha="center", va="center",
        fontsize=fontsize,
        fontfamily=FONT_FAMILY,
        fontweight="bold" if bold else "normal",
        color=fg,
        transform=ax.transAxes,
        wrap=False,
    )

# ── Header row ──────────────────────────────
header_labels = ["Parameter"] + AIRCRAFT
for c, (label, cx, cw) in enumerate(zip(header_labels, col_x, col_widths)):
    y = 1.0 - ROW_H   # top row
    draw_cell(ax, cx, y, cw, ROW_H,
              label, COL_COLORS[c], HEADER_FG,
              bold=True, fontsize=15)

# ── Data rows ───────────────────────────────
for r, param in enumerate(PARAMETERS):
    y = 1.0 - ROW_H * (r + 2)
    row_bg = ROW_COLORS[r % 2]

    # Parameter label cell
    draw_cell(ax, col_x[0], y, col_widths[0], ROW_H,
              param, row_bg, PARAM_FG, bold=True, fontsize=14)

    # Value cells
    for c, (val, cx, cw) in enumerate(zip(VALUES[r], col_x[1:], col_widths[1:])):
        draw_cell(ax, cx, y, cw, ROW_H,
                  val, row_bg, CELL_FG, bold=True, fontsize=15)

plt.tight_layout(pad=0.0)

fig.subplots_adjust(left=0, right=1, top=1, bottom=0)   # fill frame edge-to-edge

out_path = f"{OUTPUT_FILE}.{OUTPUT_FMT}"
if OUTPUT_FMT == "svg":
    matplotlib.rcParams["svg.fonttype"] = "none"   # keep text as text in SVG
    plt.savefig(out_path, format="svg", bbox_inches="tight")
else:
    plt.savefig(out_path, format="png", dpi=DPI, bbox_inches="tight",
                facecolor="white", edgecolor="none")

print(f"Saved → {out_path}")
plt.show()

"""
Redraw paper/figures/fig6_multisubband_summary.png (typeset Figure 8) from the saved result file.

The pre-registered multi-subband loop-gain probe at M <= 4, stacked for the paper's column width: (A) the truncation
sanity check (loop-gain slope s at fixed mu = 1.2 against the basis size M), (B) s(mu; M) with genuine occupancy at
the production coupling g, with the Picard boundary s = 1, (C) the pre-registered comparison statistic, the median s
over the window that survives the van Hove exclusion, against the number of occupied subbands M (and the raw median
without the exclusion). The plotted values are read, not recomputed:

  results/multisubband/r9_multisubband_sweep_data.json
      protocol_a  (rows[M, s])
      protocol_b  (g, by_M["M=1".."M=4"]: rows[mu, s], statistic[median_s_excluded, median_s_raw])

The data come from src/bdg_r9_multisubband_sweep.py. That module does not draw this figure, and no script that did is
in the repository or its history. (src/bdg_r9_multisubband_plot_M6.py extends the same three panels to M <= 6 and
writes fig7_multisubband_M1to6_extension.png.) Points of s(mu; M) with no solution (NaN) are left out of the plotted
curves, so the line joins the neighbouring points across them. This script was reconstructed from the committed
image: constants (sizes, colours, margins) were fitted until a render reproduces fig6_multisubband_summary.png pixel
for pixel.

usage (from anywhere):  py -3.13 -X utf8 -I scripts/redraw_paper_fig6.py [output.png]
default output: paper/figures/fig6_multisubband_summary.png (overwrites the committed image)
"""
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

_REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
_DATA = os.path.join(_REPO, "results", "multisubband", "r9_multisubband_sweep_data.json")
DEFAULT_OUT = os.path.join(_REPO, "paper", "figures", "fig6_multisubband_summary.png")

DPI = 140
TITLE = "The pre-registered multi-subband loop-gain probe"
FOOTER = "This is a structural result at $M\\leq4$, not a device-faithful prediction -- see text."
M_COLORS = {1: "tab:blue", 2: "tab:green", 3: "tab:red", 4: "tab:purple"}
SUPTITLE_FS = 12
TITLE_FS = 10
TOP = 0.93
HSPACE = 0.4
BOTTOM = 0.06
FOOTER_FS = 8
FOOTER_Y = 0.002
LW_A, MS_A = 1.5, 7
LW_B, MS_B, LW_PICARD, LEGEND_B_FS = 1.1, 4, 1.2, 7.5
LW_C_RED, MS_C_RED, LW_C_RAW, MS_C_RAW, LW_C_REF, LEGEND_C_FS = 1.5, 8, 1.0, 6, 1.0, 7.5


def load():
    with open(_DATA, encoding="utf-8") as f:
        return json.load(f)


def draw(data):
    pa, pb = data["protocol_a"], data["protocol_b"]

    fig, (ax, bx, cx) = plt.subplots(3, 1, figsize=(6.6, 12))

    ax.plot([r["M"] for r in pa["rows"]], [r["s"] for r in pa["rows"]], "o-", lw=LW_A, ms=MS_A)
    ax.set_xticks([r["M"] for r in pa["rows"]])
    ax.set_xlabel(r"basis size $M$ (subbands 2…$M$ unoccupied)")
    ax.set_ylabel(r"loop-gain slope $s$")
    ax.set_title("(A) Protocol (a) -- truncation sanity check:\n"
                 r"fixed $\mu=1.2$, converges smoothly (<7.6% change)", fontsize=TITLE_FS)

    for M in (1, 2, 3, 4):
        rows = pb["by_M"][f"M={M}"]["rows"]
        mu = np.array([r["mu"] for r in rows])
        s = np.array([r["s"] for r in rows])
        ok = np.isfinite(s)
        bx.plot(mu[ok], s[ok], "o-", color=M_COLORS[M], ms=MS_B, lw=LW_B, label=f"M={M}")
    bx.axhline(1.0, color="k", ls="--", lw=LW_PICARD, label=r"Picard boundary $s=1$")
    bx.set_xlabel(r"$\mu$")
    bx.set_ylabel(r"loop-gain slope $s$")
    bx.set_title("(B) Protocol (b) -- genuine occupancy, production\n"
                 fr"$g={pb['g']:.1e}$: $s=1$ crossed once $M\geq2$", fontsize=TITLE_FS)
    bx.legend(fontsize=LEGEND_B_FS, ncol=2, loc="upper right")

    Ms = [1, 2, 3, 4]
    stats = [pb["by_M"][f"M={M}"]["statistic"] for M in Ms]
    cx.plot(Ms, [st["median_s_excluded"] for st in stats], "o-", color="tab:red", ms=MS_C_RED, lw=LW_C_RED,
            label="median $s$ (van-Hove excluded, pre-registered)")
    cx.plot(Ms, [st["median_s_raw"] for st in stats], "s--", color="gray", ms=MS_C_RAW, lw=LW_C_RAW,
            label="median $s$ (raw, not excluded)")
    cx.axhline(1.0, color="k", ls=":", lw=LW_C_REF)
    cx.set_xticks(Ms)
    cx.set_xlabel(r"occupied subbands $M$")
    cx.set_ylabel(r"median loop-gain slope $s$")
    cx.set_title("(C) Pre-registered comparison statistic: sharp\n"
                 r"$M=1\to2$ jump, non-monotonic $M=2\to3\to4$", fontsize=TITLE_FS)
    cx.legend(fontsize=LEGEND_C_FS, loc="upper left")

    fig.tight_layout()
    fig.subplots_adjust(top=TOP, bottom=BOTTOM, hspace=HSPACE)
    fig.suptitle(TITLE, fontsize=SUPTITLE_FS, fontweight="bold")
    fig.text(0.5, FOOTER_Y, FOOTER, ha="center", va="bottom", fontsize=FOOTER_FS, style="italic")
    return fig


def main(out=DEFAULT_OUT):
    data = load()
    fig = draw(data)
    fig.savefig(out, dpi=DPI)
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else DEFAULT_OUT)

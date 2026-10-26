"""
Redraw results/self_consistency_bifurcation/r5_figure10_shift_vs_g.png from the saved result files.

The self-consistent boundary shift versus the feedback amplification g, in the three regimes of g: an
unambiguous boundary (g below the dead-zone onset), a dead zone at the transition, and no attracting fixed point
(above g_crit). The plotted values are read, not recomputed:

  results/self_consistency_bifurcation/r5_physical_shift_results.json
      g_sweep (g values, shifts at mu = 1.2 and 1.4), deadzone_onset_bisect, g_sweep.g_crit_est
  results/self_consistency_bifurcation/r5_results.json
      shift_comparison (the dynamically unreachable g = 3e5 root, drawn as crosses)

The data come from src/bdg_r5_physically_meaningful_shift.py (g-sweep) and src/bdg_r5_selfconsistent_ring.py
(the g = 3e5 root). Neither of them draws this figure; no script that did is in the repository or its history.
This script was reconstructed from the committed image: constants (sizes, colours, alphas, text positions) were
fitted until a render with the committed title reproduces r5_figure10_shift_vs_g.png to within 1/255 on 53 of
1.6 million pixels (all on the annotation leader line). The only change from the committed image is the title,
which no longer carries the stage prefix ("R5 Figure 10:").

usage (from anywhere):  py -3.13 -X utf8 -I scripts/redraw_r5_figure10.py [output.png]
default output: results/self_consistency_bifurcation/r5_figure10_shift_vs_g.png (overwrites the committed image)
"""
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

_REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
_DIR = os.path.join(_REPO, "results", "self_consistency_bifurcation")
DEFAULT_OUT = os.path.join(_DIR, "r5_figure10_shift_vs_g.png")

TITLE_HEAD = "Honest dose-response"
G_ROOT = 3.0e5          # coupling of the committed R5 headline root (r5_results.json, shift_comparison)
COLORS = {1.2: "#0F6E56", 1.4: "#8e44ad"}
RED, ORANGE = "#c0392b", "#e67e22"


def load():
    with open(os.path.join(_DIR, "r5_physical_shift_results.json"), encoding="utf-8") as f:
        phys = json.load(f)
    with open(os.path.join(_DIR, "r5_results.json"), encoding="utf-8") as f:
        res = json.load(f)
    return phys, res


def draw(phys, res, title_head=TITLE_HEAD):
    sweep = phys["g_sweep"]
    g = np.array(sweep["g_values"], dtype=float)
    g_crit = float(sweep["g_crit_est"])
    g_dz = float(phys["deadzone_onset_bisect"]["g_onset"])
    roots = {float(r["mu"]): float(r["selfconsistent_shift_pct"]) for r in res["shift_comparison"]}

    fig, ax = plt.subplots(figsize=(11, 6.5))
    ax.set_xscale("log")
    x_lo, x_hi = 8e2, 4e5
    ax.axvspan(x_lo, g_dz, color="#0F6E56", alpha=0.06)
    ax.axvspan(g_dz, g_crit, color=ORANGE, alpha=0.08)
    ax.axvspan(g_crit, x_hi, color=RED, alpha=0.10)

    for mu in sweep["mus"]:
        sh = np.array(sweep["shifts"][str(mu)], dtype=float)
        ok = np.isfinite(sh)
        ax.plot(g[ok], sh[ok], "o-", color=COLORS[mu], lw=2, ms=7,
                label=fr"$\mu={mu}$ (unambiguous self-consistent boundary)")
    for mu, v in roots.items():
        ax.plot([G_ROOT], [v], "X", color=COLORS[mu], ms=14, mec="k", mew=1.2, zorder=5)

    ax.axvline(g_dz, color=ORANGE, ls="--", lw=2,
               label=fr"transition dead-zone onset $\approx {g_dz:.1e}$")
    ax.axvline(g_crit, color=RED, ls="--", lw=2,
               label=fr"$g_{{\rm crit}}$ (Picard) $\approx {g_crit:.0e}$")
    ax.axhline(0, color="gray", lw=0.8)

    ax.set_xlim(x_lo, x_hi)
    ax.set_xlabel("feedback amplification $g$")
    ax.set_ylabel(r"self-consistent boundary shift $\Delta B_c/B_c$ (%)")
    ax.grid(True, which="both", alpha=0.25)

    tr = ax.get_xaxis_transform()
    for x, txt, c in ((np.sqrt(x_lo * g_dz), "regime 1\nunambiguous\nboundary", COLORS[1.2]),
                      (np.sqrt(g_dz * g_crit), "regime 2\ndead-zone at\ntransition", "#b9660f"),
                      (np.sqrt(g_crit * x_hi), "regime 3\nNO attracting\nfixed point", RED)):
        ax.text(x, 0.9558, txt, transform=tr, ha="center", va="top", color=c, fontsize=8)

    ax.annotate("g=3e5 root\n(dynamically unreachable;\nsign-flipped artifact)",
                xy=(G_ROOT, roots[1.2]), xytext=(1.1e5, 0.07123),
                color=RED, fontsize=8, arrowprops=dict(arrowstyle="-", color=RED))

    ax.legend(loc="lower left", fontsize=8.5)
    ax.set_title(title_head + r" — self-consistent shift vs coupling $g$" "\n"
                 "small, negative, grows with $g$; static answer only in regimes 1-2 "
                 "(unambiguous only in regime 1)", fontsize=11, fontweight="bold")
    # explicit axes box (pixels of the 1650 x 975 canvas): left/right/bottom are what tight_layout gives for this
    # figure; the top is fitted, because the annotation above the axes would otherwise set it differently
    fig.subplots_adjust(left=135.9791666666667 / 1650, right=1627.5 / 1650,
                        bottom=1 - 874.6833333333333 / 975, top=1 - 155.7968 / 975)
    return fig


def main(out=DEFAULT_OUT):
    phys, res = load()
    fig = draw(phys, res)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else DEFAULT_OUT)

"""
Redraw paper/figures/fig3_null_validation.png (typeset Figure 4) from the saved result file.

The validation battery for the self-consistent boundary shift in the physically meaningful regime (regime 1,
g below the dead-zone onset), stacked for the paper's column width: (a) the shift against the ring size N at g = 1e4,
which vanishes as a power law in N; (b) the shift against the bisection tolerance at g = 2e4, mu = 1.2, N = 80.
The plotted values are read, not recomputed:

  results/self_consistency_bifurcation/r5_physical_shift_results.json
      n_independence  (g, mus, Ns, shifts, power_law_exponent)
      tolerance       (g, mu, N, rows[bisect_tol, shift_pct])

The data come from src/bdg_r5_physically_meaningful_shift.py (n_independence() and tolerance_robustness()). That
module does not draw this figure, and no script that did is in the repository or its history. The dotted curves are
the power laws shift(N=40) (N/40)^p, which have the stored power_law_exponent p and pass through the first plotted
point; they are recomputed here from the stored shifts. This script was reconstructed from the committed image:
constants (sizes, colours, alphas, margins) were fitted until a render reproduces fig3_null_validation.png pixel for
pixel. The same data, drawn side by side in the style of the earlier result figures (and with the dotted curves
anchored by a least-squares prefactor instead), are in scripts/redraw_r5_figure11.py.

usage (from anywhere):  py -3.13 -X utf8 -I scripts/redraw_paper_fig3.py [output.png]
default output: paper/figures/fig3_null_validation.png (overwrites the committed image)
"""
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

_REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
_DATA = os.path.join(_REPO, "results", "self_consistency_bifurcation", "r5_physical_shift_results.json")
DEFAULT_OUT = os.path.join(_REPO, "paper", "figures", "fig3_null_validation.png")

DPI = 140
TITLE = "Validation battery for the physically meaningful coupling regime"
STYLES = {1.2: ("o-", "tab:green"), 1.4: ("D-", "tab:purple")}
TOL_COLOR = "tab:blue"
N_FIT_POINTS = 100
FIT_LW = 1.2
LEGEND_FS = 7.5
MS = 6
MS_TOL = 7
LW = 1.5
SUPTITLE_FS = 12
TITLE_FS = 10
TOP = 0.91
HSPACE = 0.35


def load():
    with open(_DATA, encoding="utf-8") as f:
        return json.load(f)


def draw(phys):
    ni, tol = phys["n_independence"], phys["tolerance"]
    Ns = np.array(ni["Ns"], dtype=float)
    N_fit = np.linspace(Ns.min(), Ns.max(), N_FIT_POINTS)

    fig, (ax, bx) = plt.subplots(2, 1, figsize=(6.6, 8.4))
    for mu in ni["mus"]:
        sh = np.array(ni["shifts"][str(mu)], dtype=float)
        p = ni["power_law_exponent"][str(mu)]
        fmt, color = STYLES[mu]
        ax.plot(Ns, sh, fmt, color=color, lw=LW, ms=MS, label=fr"$\mu={mu}$")
        ax.plot(N_fit, sh[0] * (N_fit / Ns[0]) ** p, ":", color=color, lw=FIT_LW,
                label=fr"$\mu={mu}$ fit $\propto N^{{{p:.2f}}}$")
    ax.axhline(0, color="0.6", lw=0.8)
    ax.set_xlabel("ring size $N$")
    ax.set_ylabel(r"$\Delta B_c/B_c$ (%)")
    ax.set_title(fr"(a) $N$-independence at $g={ni['g']:.0e}$ (regime 1):" "\n"
                 r"shift VANISHES as $N\to\infty$ (finite-size, not bulk)", fontsize=TITLE_FS)
    ax.legend(fontsize=LEGEND_FS, loc="lower right")

    rows = tol["rows"]
    bx.plot([r["bisect_tol"] for r in rows], [r["shift_pct"] for r in rows], "s-", color=TOL_COLOR, lw=LW, ms=MS_TOL)
    bx.set_xscale("log")
    bx.set_xlabel("bisection tolerance")
    bx.set_ylabel(r"$\Delta B_c/B_c$ (%)")
    bx.set_title(fr"(b) tolerance robustness at $g={tol['g']:.0e}$, $\mu={tol['mu']}$, $N={tol['N']}$",
                 fontsize=TITLE_FS)

    fig.tight_layout()
    fig.subplots_adjust(top=TOP, hspace=HSPACE)
    fig.suptitle(TITLE, fontsize=SUPTITLE_FS, fontweight="bold")
    return fig


def main(out=DEFAULT_OUT):
    phys = load()
    fig = draw(phys)
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    fig.savefig(out, dpi=DPI)
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else DEFAULT_OUT)

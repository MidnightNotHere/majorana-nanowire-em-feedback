"""
Redraw results/self_consistency_bifurcation/r5_figure11_headline_validation.png from the saved result file.

The validation battery for the self-consistent boundary shift in the physically meaningful regime (regime 1,
g below the dead-zone onset): (a) the shift against the ring size N at g = 1e4, which vanishes as a power law in N;
(b) the shift against the bisection tolerance at g = 2e4, mu = 1.2, N = 80. The plotted values are read, not
recomputed:

  results/self_consistency_bifurcation/r5_physical_shift_results.json
      n_independence  (g, mus, Ns, shifts, power_law_exponent)
      tolerance       (g, mu, N, rows[bisect_tol, shift_pct])

The data come from src/bdg_r5_physically_meaningful_shift.py (n_independence() and tolerance_robustness()). That
module does not draw this figure, and no script that did is in the repository or its history. The dotted curves are
the power laws -A N^p fitted to the plotted points, with p the stored power_law_exponent and A the prefactor of the
same least-squares line in log-log (np.polyfit(log N, log|shift|, 1)); they are recomputed here from the stored
shifts, and the exponents they give equal the stored ones. This script was reconstructed from the committed image:
constants (sizes, colours, alphas, margins) were fitted until a render with the committed title reproduces
r5_figure11_headline_validation.png pixel for pixel. The only change from the committed image is the title, which
no longer carries the stage prefix ("R5 Figure 11:").

usage (from anywhere):  py -3.13 -X utf8 -I scripts/redraw_r5_figure11.py [output.png]
default output: results/self_consistency_bifurcation/r5_figure11_headline_validation.png (overwrites the committed
image)
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
DEFAULT_OUT = os.path.join(_DIR, "r5_figure11_headline_validation.png")

TITLE_HEAD = "Validation battery in the physically meaningful regime"
COLORS = {1.2: "#0F6E56", 1.4: "#8e44ad"}
TOL_COLOR = "#2471a3"
N_FIT_POINTS = 100
FIT_LW = 1.4
LEGEND_FS = 9


def _fit_grid(Ns):
    return np.linspace(Ns.min(), Ns.max(), N_FIT_POINTS)


def load():
    with open(os.path.join(_DIR, "r5_physical_shift_results.json"), encoding="utf-8") as f:
        return json.load(f)


def draw(phys, title_head=TITLE_HEAD):
    ni, tol = phys["n_independence"], phys["tolerance"]
    Ns = np.array(ni["Ns"], dtype=float)
    N_fit = _fit_grid(Ns)

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(14, 5.5))
    for mu in ni["mus"]:
        sh = np.array(ni["shifts"][str(mu)], dtype=float)
        p = ni["power_law_exponent"][str(mu)]
        slope, icpt = np.polyfit(np.log(Ns), np.log(-sh), 1)
        ax.plot(Ns, sh, "o-", color=COLORS[mu], lw=2, ms=8, label=fr"$\mu={mu}$")
        ax.plot(N_fit, -np.exp(icpt) * N_fit ** slope, ":", color=COLORS[mu], lw=FIT_LW,
                label=fr"$\mu={mu}$ fit $\propto N^{{{p:.2f}}}$")
    ax.axhline(0, color="gray", lw=0.8)
    ax.set_xlabel("ring size $N$")
    ax.set_ylabel(r"$\Delta B_c/B_c$ (%)")
    ax.set_title(fr"(a) N-independence at $g={ni['g']:.0e}$ (regime 1):" "\n"
                 r"shift VANISHES as $N\to\infty$ (finite-size, not bulk)", fontsize=11)
    ax.legend(fontsize=LEGEND_FS, loc="lower right")
    ax.grid(True, alpha=0.3)

    rows = tol["rows"]
    bx.plot([r["bisect_tol"] for r in rows], [r["shift_pct"] for r in rows], "s-", color=TOL_COLOR, lw=2, ms=8)
    bx.set_xscale("log")
    bx.invert_xaxis()
    bx.set_xlabel("bisection tolerance")
    bx.set_ylabel(r"$\Delta B_c/B_c$ (%)")
    bx.set_title(fr"(b) tolerance robustness at $g={tol['g']:.0e}$, $\mu={tol['mu']}$, $N={tol['N']}$", fontsize=11)
    bx.grid(True, alpha=0.3)

    fig.suptitle(title_head + r" (regime 1, $g<g_{\rm deadzone}$)", fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    return fig


def main(out=DEFAULT_OUT):
    phys = load()
    fig = draw(phys)
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else DEFAULT_OUT)

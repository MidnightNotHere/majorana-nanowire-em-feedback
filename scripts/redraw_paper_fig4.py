"""
Redraw paper/figures/fig4_bounded_pole.png from the saved reactive-impedance pole data.

Nothing is recomputed: the plotted curves are read from
results/reactive_impedance/pole_characterization/r9_pole_data.json. The two panels are the
Figure-1 panels of bdg_r9_pole_characterization.make_figures(), stacked vertically at
single-column width under a shorter super-title so that the title is not clipped when the
figure is scaled to the journal column. The plotting calls, constants and data are those of
make_figures(); only the figure geometry and the super-title differ.

usage:  py -3.13 scripts/redraw_paper_fig4.py [output.png]
default output: paper/figures/fig4_bounded_pole.png
"""
import json
import os
import sys

_REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(_REPO, "src"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import bdg_r9_pole_characterization as pole

# Importing the module applies its own rcParams (large fonts, thick lines); the paper
# figure uses matplotlib's defaults.
matplotlib.rcdefaults()

DATA_FILE = pole._out("r9_pole_data.json")
DEFAULT_OUT = os.path.join(_REPO, "paper", "figures", "fig4_bounded_pole.png")


def draw(data, out):
    br = [r for r in data["branch"] if np.isfinite(r["s"]) and r["s"] < 1.0]
    g = np.array([r["g"] for r in br]); enh = np.array([r["enh"] for r in br])
    det = data["detach"]; fit = data["fit"]; gcr = data["gcrit"]

    fig, (axA, axB) = plt.subplots(2, 1, figsize=(6.6, 8.6))
    axA.plot(g, enh, "o-", color="#2ca02c", ms=4, label=r"$1/(1-s)$ principal branch")
    if np.isfinite(det.get("g_detach", np.nan)):
        axA.axvline(det["g_detach"], color="#d62728", ls="--", lw=1.3,
                    label=f"detachment g≈{det['g_detach']:.2e}\n(enh caps at {det['enh_at_detach']:.1f})")
    axA.axvspan(pole.G_CRIT_R5_LO, pole.G_CRIT_R5_HI, color="orange", alpha=0.18,
                label="g_crit band (1.1–1.2e5)")
    axA.axvline(pole.G_R7, color="k", ls=":", lw=1.1, label=f"production g={pole.G_R7:.1e}")
    if gcr.get("g_crit_slope_cross"):
        axA.axvline(gcr["g_crit_slope_cross"], color="purple", ls="-.", lw=1.1,
                    label=f"slope=1 rule cross ≈{gcr['g_crit_slope_cross']:.2e}")
    axA.set_xlabel("feedback strength  g"); axA.set_ylabel(r"enhancement $1/(1-s)$")
    axA.set_title("(A)  Enhancement caps at ~%.1f (N=80), then detaches\n"
                  "(fluxoid sawtooth) — s never reaches 1"
                  % fit.get("max_enh_attained", float("nan")), fontsize=10)
    axA.legend(fontsize=7.5, loc="upper left")

    la = sorted(data["nscale"]["L_axis"], key=lambda r: r["N"])
    fl = sorted(data["nscale"]["fixedL"], key=lambda r: r["N"])
    axB.plot([r["N"] for r in la], [r["max_enh"] for r in la], "o-",
             color="#1f77b4", label="L-axis (a=10 nm; L & N grow)")
    axB.plot([r["N"] for r in fl], [r["max_enh"] for r in fl], "s-",
             color="#d62728", label=f"fixed-L mesh (L={data['nscale']['fixedL_um']} μm)")
    axB.axhline(1.0, color="0.7", lw=0.8)
    axB.set_xlabel("N (sites)"); axB.set_ylabel("max enhancement reached")
    axB.set_title("(B)  Cap is N-dependent: it CLIMBS as the mesh\n"
                  "refines (pole sharpens) — bound is not universal", fontsize=10)
    axB.legend(fontsize=8)

    fig.suptitle("The reactive-impedance 'pole' near $g_{\\rm crit}$", fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.subplots_adjust(top=0.91, hspace=0.35)
    fig.savefig(out, dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_OUT
    with open(DATA_FILE, encoding="utf-8") as f:
        data = pole._coerce(json.load(f))
    draw(data, out)
    print("Saved:", os.path.abspath(out))

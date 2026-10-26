"""
Redraw results/multisubband/r9_multisubband_summary.png from the saved result file.

The pre-registered multi-subband loop-gain probe at M <= 4, in three panels side by side: (A) the truncation sanity
check (loop-gain slope s at fixed mu = 1.2 against the basis size M, with the single-channel baseline s ~ 0.29),
(B) s(mu; M) with genuine occupancy at the production coupling g, with the Picard boundary s = 1 and a cross at the
bottom for every mu with no root (detached), (C) the pre-registered statistic, the median s over the window that
survives the van Hove exclusion, against the number of occupied subbands M (and the raw median without the exclusion),
each point annotated with its surviving / finite / pre-registered point counts. The plotted values are read, not
recomputed:

  results/multisubband/r9_multisubband_sweep_data.json
      protocol_a  (rows[M, s])
      protocol_b  (g, by_M["M=1".."M=4"]: rows[mu, s], statistic[median_s_excluded, median_s_raw,
                   n_surviving, n_raw_finite])

The data come from src/bdg_r9_multisubband_sweep.py. That module does not draw this figure, and no script that did is
in the repository or its history. This script was reconstructed from the committed image: constants (sizes, colours,
margins, text positions) were fitted until a render with the committed text reproduces r9_multisubband_summary.png
to within one grey level (panels A and B are identical; in panel C 2,952 of the 1,834,560 pixels, around the markers,
differ by 1/255). The baseline level s = 0.29 of panel A is a drawing constant, as in the committed image. The text
that the committed image carried differently, and that is changed here, is the stage
prefix of the figure title ("R9 Track 4 --"), the stage code of the baseline legend entry ("R9-T1 single-channel
baseline"), and the methods-note file name, now given with its path (docs/methods/r9_multisubband_methods.md).
paper/figures/fig6_multisubband_summary.png (typeset Figure 8) is the same figure laid out in a single column; it is
redrawn by scripts/redraw_paper_fig6.py.

usage (from anywhere):  py -3.13 -X utf8 -I scripts/redraw_r9_multisubband_summary.py [output.png]
default output: results/multisubband/r9_multisubband_summary.png (overwrites the committed image)
"""
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

_REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
_DIR = os.path.join(_REPO, "results", "multisubband")
DEFAULT_OUT = os.path.join(_DIR, "r9_multisubband_summary.png")

DPI = 140
SUPTITLE_HEAD = "Multi-subband loop-gain probe"
NOTE_REF = "docs/methods/r9_multisubband_methods.md"
BASELINE_LABEL = "single-channel baseline (s≈0.29)"
BASELINE_S = 0.29
M_COLORS = {1: "tab:blue", 2: "tab:green", 3: "tab:red", 4: "tab:purple"}
SUPTITLE_FS = 11
TITLE_FS = 10
LABEL_FS = 10
TIGHT_RECT = [0, 0, 1, 0.9]
LW_A, MS_A, LEGEND_A_FS, LW_BASE, BASE_ALPHA = 1.5, 8, 8, 1.2, 1.0
LW_B, MS_B, LW_S1, COLOR_S1, LEGEND_B_FS, X_Y, X_MS = 1.5, 5, 1.2, "0.4", 8, -0.3, 6
LW_C_RED, MS_C_RED, LW_C_RAW, MS_C_RAW, LEGEND_C_FS, NOTE_FS = 2, 9, 1.5, 7, 8, 7.5
NOTE_OFFSET = (6, 8)
LW_C_REF, COLOR_C_REF = 1.2, "0.4"


def load():
    with open(os.path.join(_DIR, "r9_multisubband_sweep_data.json"), encoding="utf-8") as f:
        return json.load(f)


def draw(data, head=SUPTITLE_HEAD, note_ref=NOTE_REF, baseline_label=BASELINE_LABEL):
    pa, pb = data["protocol_a"], data["protocol_b"]

    fig, (ax, bx, cx) = plt.subplots(1, 3, figsize=(18, 5.2))

    ax.plot([r["M"] for r in pa["rows"]], [r["s"] for r in pa["rows"]], "o-", lw=LW_A, ms=MS_A)
    ax.axhline(BASELINE_S, color="gray", ls=":", lw=LW_BASE, alpha=BASE_ALPHA, label=baseline_label)
    ax.set_xticks([r["M"] for r in pa["rows"]])
    ax.set_xlabel("M (basis size, subbands 2..M unoccupied)", fontsize=LABEL_FS)
    ax.set_ylabel("s", fontsize=LABEL_FS)
    ax.set_title("(A) Protocol (a): truncation sanity check\n"
                 "fixed μ=1.2 — converges smoothly (<7.6% change)", fontsize=TITLE_FS)
    ax.legend(fontsize=LEGEND_A_FS, loc="center right")

    for M in (1, 2, 3, 4):
        rows = pb["by_M"][f"M={M}"]["rows"]
        mu = np.array([r["mu"] for r in rows])
        s = np.array([r["s"] for r in rows])
        ok = np.isfinite(s)
        bx.plot(mu[ok], s[ok], "o-", color=M_COLORS[M], lw=LW_B, ms=MS_B, label=f"M={M} (n={len(rows)})")
        bx.plot(mu[~ok], np.full((~ok).sum(), X_Y), "x", color=M_COLORS[M], ms=X_MS)
    bx.axhline(1.0, color=COLOR_S1, ls="--", lw=LW_S1, label="s=1 (Picard-stability boundary)")
    bx.set_xlabel(r"$\mu$ (Δ-units)", fontsize=LABEL_FS)
    bx.set_ylabel("loop-gain slope $s$", fontsize=LABEL_FS)
    bx.set_title("(B) Protocol (b): genuine occupancy, s(μ;M)\n"
                 f"at production g={pb['g'] / 1e4:g}e4 (×'s at bottom = detached/no root)", fontsize=TITLE_FS)
    bx.legend(fontsize=LEGEND_B_FS, loc="upper right")

    Ms = [1, 2, 3, 4]
    stats = [pb["by_M"][f"M={M}"]["statistic"] for M in Ms]
    med_ex = [st["median_s_excluded"] for st in stats]
    med_raw = [st["median_s_raw"] for st in stats]
    cx.plot(Ms, med_ex, "o-", color="tab:red", lw=LW_C_RED, ms=MS_C_RED,
            label="median (van Hove excluded) — PRE-REGISTERED")
    cx.plot(Ms, med_raw, "s--", color="gray", lw=LW_C_RAW, ms=MS_C_RAW, label="median (raw, non-excluded)")
    cx.axhline(1.0, color=COLOR_C_REF, ls=":", lw=LW_C_REF)
    for M, st, y in zip(Ms, stats, med_ex):
        n_total = len(pb["by_M"][f"M={M}"]["rows"])
        cx.annotate(f"n={st['n_surviving']}/{st['n_raw_finite']} surv.\n(of {n_total} pre-reg.)", (M, y),
                    textcoords="offset points", xytext=NOTE_OFFSET, fontsize=NOTE_FS)
    cx.set_xticks(Ms)
    cx.set_xlabel("M (occupied subbands)", fontsize=LABEL_FS)
    cx.set_ylabel("median $s$ over surviving window", fontsize=LABEL_FS)
    cx.set_title("(C) The pre-registered statistic:\n"
                 "sharp M=1→2 jump, then non-monotonic M=2..4\n"
                 "(M=4 densified 5→15 pts — see §8.2b)", fontsize=TITLE_FS)
    cx.legend(fontsize=LEGEND_C_FS, loc="upper left")

    fig.suptitle(f"{head}: s crosses the Picard-stability boundary at production g once M≥2\n"
                 f"(STRUCTURAL result, M≤4 — NOT a device-faithful prediction; see {note_ref} §8-9)",
                 fontsize=SUPTITLE_FS, fontweight="bold")
    fig.tight_layout(rect=TIGHT_RECT)
    return fig


def main(out=DEFAULT_OUT):
    data = load()
    fig = draw(data)
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    fig.savefig(out, dpi=DPI)
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else DEFAULT_OUT)

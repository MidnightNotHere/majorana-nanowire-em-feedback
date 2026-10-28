"""
R9 (Track 4 extension): figure for the M=1..6 loop-gain probe.
=============================================================
Combines the original M=1..4 pre-registered sweep
(`results/multisubband/r9_multisubband_sweep_data.json`) with the M=5,6
extension (`r9_multisubband_sweep_M5_M6.json`, `r9_multisubband_protocol_a_M6.json`)
into a single 3-panel figure, in the same layout as the original
fig6_multisubband_summary.png (A: truncation convergence, B: s(mu;M) curves,
C: pre-registered median statistic vs M). This is a NEW figure file --
the original fig6 is left untouched, so the M<=4 result and its provenance
remain independently reproducible.
"""
import json
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# repo root, resolved relative to THIS file so the script runs from any checkout
BASE = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
R9 = os.path.join(BASE, "results", "multisubband")

orig = json.load(open(os.path.join(R9, "r9_multisubband_sweep_data.json")))
ext = json.load(open(os.path.join(R9, "r9_multisubband_sweep_M5_M6.json")))
pa6 = json.load(open(os.path.join(R9, "r9_multisubband_protocol_a_M6.json")))

fig, axes = plt.subplots(3, 1, figsize=(4.4, 10.8))

# --- Panel A: truncation convergence (fixed occupancy sanity check) ---
ax = axes[0]
rows_a = pa6["rows"]
Ms_a = [r["M"] for r in rows_a]
s_a = [r["s"] for r in rows_a]
ax.plot(Ms_a, s_a, "o-", color="tab:blue")
ax.set_xlabel("M (basis size, fixed occupancy)")
ax.set_ylabel(r"$s$")
ax.set_title("(A) Truncation convergence\n(sanity check, not physics)")
ax.set_xticks(Ms_a)
ax.grid(alpha=0.3)

# --- Panel B: s(mu;M) curves, genuine occupancy ---
ax = axes[1]
colors = plt.cm.viridis(np.linspace(0, 0.9, 6))
for i, M in enumerate(range(1, 7)):
    if M <= 4:
        rows = orig["protocol_b"]["by_M"][f"M={M}"]["rows"]
    else:
        rows = ext["sweep"]["by_M"][f"M={M}"]["rows"]
    rows = sorted(rows, key=lambda r: r["mu"])
    mus = np.array([r["mu"] for r in rows])
    svals = np.array([r["s"] for r in rows])
    # keep NaNs in place (matplotlib breaks the line there) rather than
    # filtering them out, which would draw a misleading direct line straight
    # through a genuine no-solution (detachment) gap
    ax.plot(mus, svals, "o-", color=colors[i], label=f"M={M}", ms=3)
ax.axhline(1.0, color="k", ls="--", lw=1, label=r"$s=1$")
ax.set_xlabel(r"$\mu$ ($\Delta$-units)")
ax.set_ylabel(r"$s(\mu)$")
ax.set_yscale("symlog", linthresh=1e-3)
ax.set_title("(B) Genuine-occupancy sweep\n$g=2.9\\times10^4$")
ax.legend(fontsize=6, ncol=2, loc="upper right")
ax.grid(alpha=0.3)

# --- Panel C: pre-registered median statistic vs M ---
ax = axes[2]
Ms_c = list(range(1, 7))
med = []
detach_rate = []
for M in Ms_c:
    if M <= 4:
        st = orig["protocol_b"]["by_M"][f"M={M}"]["statistic"]
        n_total = len(orig["protocol_b"]["by_M"][f"M={M}"]["rows"])
    else:
        st = ext["sweep"]["by_M"][f"M={M}"]["statistic"]
        n_total = len(ext["sweep"]["by_M"][f"M={M}"]["rows"])
    med.append(st["median_s_excluded"])
    n_detach = n_total - st["n_raw_finite"]
    detach_rate.append(100.0 * n_detach / n_total)

ax.plot(Ms_c, med, "o-", color="tab:red", label="median $s$ (excluded)")
ax.axhline(1.0, color="k", ls="--", lw=1)
ax.set_xlabel("M (occupied subbands)")
ax.set_ylabel("median $s$", color="tab:red")
ax.set_yscale("symlog", linthresh=1e-4)
ax.set_xticks(Ms_c)
ax.tick_params(axis="y", labelcolor="tab:red")
ax.set_title("(C) Pre-registered statistic\nvs $M$")
ax.grid(alpha=0.3)

ax2 = ax.twinx()
ax2.plot(Ms_c, detach_rate, "s--", color="tab:gray", alpha=0.7, label="detachment rate")
ax2.set_ylabel("detachment rate (%)", color="tab:gray")
ax2.tick_params(axis="y", labelcolor="tab:gray")
ax2.set_ylim(0, 60)

plt.tight_layout()
outpath = os.path.join(R9, "fig7_multisubband_M1to6_extension.png")
os.makedirs(R9, exist_ok=True)
plt.savefig(outpath, dpi=200)
print("wrote", outpath)

# also copy into the paper's figures/ directory for \includegraphics
import shutil
paper_figures = os.path.join(BASE, "paper", "figures")
os.makedirs(paper_figures, exist_ok=True)
shutil.copy(outpath, os.path.join(paper_figures, "fig7_multisubband_M1to6_extension.png"))
print("copied into paper/figures/")

# print the exact numbers used, for the caption/text
print("\nM : median_s_excluded : n_surviving/n_raw_finite : detachment_rate%")
for M, m, r in zip(Ms_c, med, detach_rate):
    print(f"{M} : {m:.6g} : {r:.1f}%")

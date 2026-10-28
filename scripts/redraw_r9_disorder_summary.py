"""
Redraw results/disorder/r9_disorder_summary.png from the saved result file, without re-running the disorder ensembles.

The figure is drawn by make_figures() in src/bdg_r9_disorder.py, which the full run (bdg_r9_disorder.py --full, 4,627 s)
calls after the ensembles have finished. This script calls the same function on the saved data, so the
plotted values are read, not recomputed:

  results/disorder/r9_disorder_data.json
      task4_1_boundary_smearing, task4_2_subgap_states, task4_3_null_vs_W  (rows per disorder strength W)

The image it writes is identical, pixel for pixel, to the committed results/disorder/r9_disorder_summary.png, and to
paper/figures/fig5_disorder_summary.png (typeset Figure 6), which is a copy of it. It does not redraw
r9_disorder_summary.SUPERSEDED_task4_3.png: that image was drawn from r9_disorder_data.SUPERSEDED_task4_3.json by an
earlier version of make_figures() and is kept only as the record of the superseded result.

usage (from anywhere):  py -3.13 -X utf8 -I scripts/redraw_r9_disorder_summary.py [output.png]
default output: results/disorder/r9_disorder_summary.png (overwrites the committed image)
"""
import json
import os
import shutil
import sys
import tempfile

_REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.join(_REPO, "src"))

import bdg_r9_disorder  # noqa: E402

_DATA = os.path.join(_REPO, "results", "disorder", "r9_disorder_data.json")
DEFAULT_OUT = os.path.join(_REPO, "results", "disorder", "r9_disorder_summary.png")


def main(out=DEFAULT_OUT):
    with open(_DATA, encoding="utf-8") as f:
        data = json.load(f)
    # make_figures() writes into its own results folder; draw into a scratch folder and move the image to `out`
    with tempfile.TemporaryDirectory() as tmp:
        bdg_r9_disorder._out = lambda name: os.path.join(tmp, os.path.basename(str(name)))
        (drawn,) = bdg_r9_disorder.make_figures(data)
        os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
        shutil.copyfile(drawn, out)
    print("wrote", out)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else DEFAULT_OUT)

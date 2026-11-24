"""BoxGraph software adapter invoking getlinefracinsidedroplet() directly from VASPDroplet_simulate.ipynb.

Demonstrates:
  - scenario='edfig6_11seg': 12 beads (11 segments of 100 nm = 1.1 um) spanning a 1.0 um droplet with 100 nm
    protruding on one side -> returns 10/11 = 0.90909 (matching Extended Data Fig. 6A's ~0.9092).
  - scenario='edfig6_12seg': 13 beads (12 segments of 100 nm = 1.2 um) centered in a 1.0 um droplet with 100 nm
    protruding at each end -> returns 10/12 = 0.83333.
  - scenario='single_chord_two_intersect': 2 beads (1 segment of 1.2 um crossing a 1.0 um droplet twice) ->
    exposes lines 313-314 bug (overwriting leninside with = and swapping inside/outside fractions),
    returning 200/1200 = 0.16667 instead of the true wetted fraction 1000/1200 = 0.83333.
"""

import json
import os
from pathlib import Path
import numpy as np

_NS = None


def _load_notebook():
    global _NS
    if _NS is not None:
        return _NS
    nb_path = Path(os.environ["BOXGRAPH_SOURCE"]) / "VASPDroplet_simulate.ipynb"
    nb = json.loads(nb_path.read_text(encoding="utf-8"))
    ns = {}
    for idx in (1, 3, 5, 7):
        src = "".join(nb["cells"][idx]["source"])
        if idx == 1:
            src = "\n".join(ln for ln in src.splitlines() if "matplotlib" not in ln)
        exec(compile(src, f"VASPDroplet_simulate.ipynb#cell{idx}", "exec"), ns)
    _NS = ns
    return _NS


def run(inputs, settings):
    ns = _load_notebook()
    sp = ns["SystemParams"]()
    sp.dimval = 2
    r_nm = float(inputs["D_droplet"]) * 500.0
    lcyl = float(inputs["lcyl"])
    com = float(inputs["ellipse_com"])
    sp.ellipse_com = [com, com]
    sp.Vdropfactor = r_nm * r_nm
    sp.lcyl = lcyl

    scenario = str(settings.get("scenario", "edfig6_11seg"))
    if scenario == "edfig6_11seg":
        # 12 beads = 11 segments of 100 nm = 1.1 um: from -R to +R + 100 nm
        n_seg = int(round((2.0 * r_nm + lcyl) / lcyl))
        sp.Nbeads = n_seg + 1
        xs = np.linspace(com - r_nm, com + r_nm + lcyl, sp.Nbeads)
    elif scenario == "edfig6_12seg":
        # 13 beads = 12 segments of 100 nm = 1.2 um: from -R - 100 nm to +R + 100 nm
        n_seg = int(round((2.0 * r_nm + 2.0 * lcyl) / lcyl))
        sp.Nbeads = n_seg + 1
        xs = np.linspace(com - r_nm - lcyl, com + r_nm + lcyl, sp.Nbeads)
    elif scenario == "single_chord_two_intersect":
        # 2 beads = 1 straight segment of length 2*R + 200 nm intersecting circle twice
        sp.Nbeads = 2
        xs = np.array([com - r_nm - lcyl, com + r_nm + lcyl])
    else:
        raise ValueError(f"Unknown scenario: {scenario!r}")

    coords = np.column_stack([xs, np.full(sp.Nbeads, com)])
    opt_var = np.concatenate([[r_nm], coords.ravel()])
    frac = float(ns["getlinefracinsidedroplet"](opt_var, sp))
    return {"wetted_fraction": frac}

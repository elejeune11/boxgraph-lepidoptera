"""BoxGraph software adapter invoking getsurfaceF() directly from VASPDroplet_simulate.ipynb."""

import json
import math
import os
from pathlib import Path

_NS = None


def _load_notebook():
    global _NS
    if _NS is not None:
        return _NS
    nb_path = Path(os.environ["BOXGRAPH_SOURCE"]) / "VASPDroplet_simulate.ipynb"
    nb = json.loads(nb_path.read_text(encoding="utf-8"))
    ns = {}
    for idx in (1, 3, 5, 7, 9):
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
    r_nm = float(inputs["R"])
    ar = float(inputs["AR"])
    sp.Vdropfactor = r_nm * r_nm
    sp.surface_T = [float(inputs["gamma_I"]), 0.0, 0.0]
    a = r_nm * math.sqrt(ar)
    val = ns["getsurfaceF"]([a], sp, False)
    return {"E_surf": float(val)}

"""BoxGraph software adapter invoking getsurfaceTension() directly from VASPDroplet_simulate.ipynb."""

import json
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
    sp.surface_T = [0.0, float(inputs["gamma_AVR"]), float(inputs["gamma_AVD"])]
    val = ns["getsurfaceTension"](float(inputs["e"]), sp)
    return {"gamma": float(val)}

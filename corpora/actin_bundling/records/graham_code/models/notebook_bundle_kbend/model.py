"""BoxGraph software adapter invoking getkbend() directly from VASPDroplet_simulate.ipynb."""

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
    sp.Nfil = float(inputs["Nfil"])
    sp.lcyl = float(inputs["lcyl"])
    sp.kflex_rigidity_f = float(inputs["kflex_fil"])
    sp.alpha = float(inputs["alpha"])
    sp.nonlinearbool = bool(settings.get("nonlinearbool", False))
    kbend = float(ns["getkbend"](sp))
    return {"k_bend": kbend, "k_flex_bundle": kbend * sp.lcyl}

"""Live software adapter executing the Arp2/3 molecule count line from Cell 15 of CytoSimanalyze_VASP_Arp_VCA.ipynb."""

import json
import os
from pathlib import Path
import numpy as np
import scipy.constants

_EXPR = None


def _load_expr():
    global _EXPR
    if _EXPR is not None:
        return _EXPR
    nb_path = Path(os.environ["BOXGRAPH_SOURCE"]) / "CytoSimanalyze_VASP_Arp_VCA.ipynb"
    nb = json.loads(nb_path.read_text(encoding="utf-8"))
    for ln in nb["cells"][15]["source"]:
        s = ln.strip()
        if s.startswith("Narp ="):
            _EXPR = compile(s, "CytoSimanalyze_VASP_Arp_VCA.ipynb#cell15", "exec")
            return _EXPR
    raise RuntimeError("Narp assignment not found in Cell 15")


def run(inputs, settings):
    code = _load_expr()
    ns = {"np": np, "scipy": scipy, "CArp": np.array([float(inputs["C_Arp"])])}
    exec(code, ns)
    return {"N_arp": float(ns["Narp"][0])}

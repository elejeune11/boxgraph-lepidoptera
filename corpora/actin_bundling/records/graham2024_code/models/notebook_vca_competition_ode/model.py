"""Live software adapter executing model_VCA() directly from Cell 23 of CytoSimanalyze_VASP_Arp_VCA.ipynb."""

import ast
import json
import math
import os
from pathlib import Path
import numpy as np
import scipy.constants
from scipy.integrate import odeint

_MODEL_FN = None


def _load_model_fn():
    global _MODEL_FN
    if _MODEL_FN is not None:
        return _MODEL_FN
    nb_path = Path(os.environ["BOXGRAPH_SOURCE"]) / "CytoSimanalyze_VASP_Arp_VCA.ipynb"
    nb = json.loads(nb_path.read_text(encoding="utf-8"))
    src = "".join(nb["cells"][23]["source"])
    tree = ast.parse(src)
    fn_node = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "model_VCA")
    mod = ast.Module(body=[fn_node], type_ignores=[])
    ns = {"np": np}
    exec(compile(mod, "CytoSimanalyze_VASP_Arp_VCA.ipynb#cell23", "exec"), ns)
    _MODEL_FN = ns["model_VCA"]
    return _MODEL_FN


def run(inputs, settings):
    model_fn = _load_model_fn()
    Nfil = float(inputs["N_fil"])
    TFactin = float(inputs["TFactin"])
    grow = [float(inputs["k_grow_plus"]), 0.0]
    shrink = [0.0, 0.0]
    vcakin = [float(inputs["k_bind_vca"]), float(inputs["k_unbind_vca"])]
    nvca = float(inputs["N_vca"])
    L0 = float(inputs["L0_per_fil"]) * Nfil
    R = float(inputs["R_drop"])
    l_mono = float(inputs["l_mono"])

    t_max = float(settings.get("t_grid_max", 600.0))
    n_grid = int(settings.get("n_grid", 601))
    eval_idx = int(settings.get("eval_index", 600))
    rtol = float(settings.get("rtol", 1.0e-8))
    atol = float(settings.get("atol", 1.0e-8))

    t = np.linspace(0.0, t_max, n_grid)
    y0 = odeint(model_fn, np.array([L0, 0.0, 0.0, 0.0, 0.0]), t, args=(Nfil, TFactin, grow, shrink, vcakin), rtol=rtol, atol=atol)
    y = odeint(model_fn, np.array([L0, nvca, 0.0, 0.0, 0.0]), t, args=(Nfil, TFactin, grow, shrink, vcakin), rtol=rtol, atol=atol)

    L_no_vca = float(y0[eval_idx, 0])
    L_val = float(y[eval_idx, 0])
    L_vca = float(y[eval_idx, 4])

    vol_L = (4.0 / 3.0) * math.pi * (R ** 3) * 1e-15
    uM_per_um = (1.0 / l_mono) / (scipy.constants.Avogadro * vol_L) * 1e6

    return {
        "L_total": L_val,
        "L_vca": L_vca,
        "F_actin_conc": L_val * uM_per_um,
        "delta_F_actin": (L_no_vca - L_val) * uM_per_um,
        "frac_F_actin_reduction": 1.0 - L_val / L_no_vca,
        "frac_F_actin_remaining": L_val / L_no_vca,
    }

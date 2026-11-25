"""Live software adapter executing model() directly from Cell 21 of CytoSimanalyze_VASP_Arp_VCA.ipynb."""

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
    src = "".join(nb["cells"][21]["source"])
    tree = ast.parse(src)
    fn_node = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "model")
    mod = ast.Module(body=[fn_node], type_ignores=[])
    ns = {"np": np}
    exec(compile(mod, "CytoSimanalyze_VASP_Arp_VCA.ipynb#cell21", "exec"), ns)
    _MODEL_FN = ns["model"]
    return _MODEL_FN


def run(inputs, settings):
    model_fn = _load_model_fn()
    Nfil = float(inputs["N_fil"])
    TFactin = float(inputs["TFactin"])
    grow = [float(inputs["k_grow_plus"]), float(settings.get("k_grow_minus", 0.0))]
    shrink = [float(settings.get("k_shrink_plus", 0.0)), float(settings.get("k_shrink_minus", 0.0))]
    L0 = float(inputs["L0_per_fil"]) * Nfil
    R = float(inputs["R_drop"])
    l_mono = float(inputs["l_mono"])

    t_max = float(settings.get("t_grid_max", 600.0))
    n_grid = int(settings.get("n_grid", 601))
    eval_idx = int(settings.get("eval_index", 600))
    rtol = float(settings.get("rtol", 1.0e-8))
    atol = float(settings.get("atol", 1.0e-8))

    t = np.linspace(0.0, t_max, n_grid)
    y = odeint(model_fn, np.array([L0, 0.0, 0.0]), t, args=(Nfil, TFactin, grow, shrink), rtol=rtol, atol=atol)

    dt = np.diff(t)
    dPdt = np.diff(y[:, 1]) / dt
    L_val = float(y[eval_idx, 0])
    dP_val = float(dPdt[eval_idx]) if eval_idx < len(dPdt) else float(dPdt[-1])
    tm_idx = 1000 if len(dPdt) > 1000 else len(dPdt) - 1
    tm_rate = float(dPdt[tm_idx] / Nfil)

    vol_L = (4.0 / 3.0) * math.pi * (R ** 3) * 1e-15
    uM_per_um = (1.0 / l_mono) / (scipy.constants.Avogadro * vol_L) * 1e6

    return {
        "L_total": L_val,
        "F_actin_conc": L_val * uM_per_um,
        "dP_dt": dP_val,
        "treadmilling_rate": tm_rate,
    }

import importlib.util
import math
from pathlib import Path
import numpy as np

_SSC_PATH = Path(__file__).resolve().parent.parent / "steady_state_continuation_eq7" / "model.py"
_SPEC = importlib.util.spec_from_file_location("_bg_steady_state_continuation_eq7", _SSC_PATH)
_SSC_MOD = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_SSC_MOD)
continue_steady = _SSC_MOD.continue_steady


def run(inputs, settings):
    x = float(inputs["x"])
    a = inputs["a"]
    a_Delta = inputs["a_Delta"]
    kappa = inputs["kappa"]
    Q = inputs["Q"]
    rho = inputs["rho"]
    alpha = inputs["alpha"]
    ell = inputs["ell"]
    p1_bar = inputs["phi1_bar"]
    p2_bar = inputs["phi2_bar"]
    n_per = max(1, int(round(inputs["n_period"])))
    amp1 = inputs["init_amp_1"]
    amp2 = inputs["init_amp_2"]

    M = int(settings.get("n_modes", 36))
    tol = float(settings.get("tol", 1e-11))
    x0 = float(settings.get("phase_x0", 0.0))

    u, ok = continue_steady(
        a, a_Delta, kappa, Q, rho, alpha, ell, p1_bar, p2_bar, n_per, M, tol, amp1, amp2, settings
    )
    if not ok:
        return {"phi1": None, "phi2": None}

    c1, c2 = u[:M], u[M:]
    ms = np.arange(1, M + 1)
    theta = 2.0 * math.pi * n_per * (x - x0)
    cos_v = np.cos(theta * ms)
    return {
        "phi1": float(p1_bar + np.dot(cos_v, c1)),
        "phi2": float(p2_bar + np.dot(cos_v, c2)),
    }

import importlib.util
from pathlib import Path
import numpy as np

_base_path = Path(__file__).resolve().parent.parent / "sphere_hydrodynamics_single_step_eqs32_33" / "model.py"
_spec = importlib.util.spec_from_file_location("_sphere_hydro_base_s2a", _base_path)
_base = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_base)


def compute_orders(inputs, settings, pole_kappa_coeff=8.0):
    N_base = int(round(settings.get("N_base", 51)))
    if (N_base - 1) % 2 != 0:
        N_base += 1
    N_ref = int(round(settings.get("N_ref", 801)))
    # Ensure nested grids N_base, 2*N_base-1, 4*N_base-3 divide N_ref-1
    n1 = N_base
    n2 = 2 * (n1 - 1) + 1
    n3 = 4 * (n1 - 1) + 1
    mult = max(2, int(round((N_ref - 1) / (n3 - 1))))
    N_ref_exact = mult * (n3 - 1) + 1

    vs_ref, vn_ref, _ = _base.solve_single_step(inputs, {"N": N_ref_exact}, pole_kappa_coeff=pole_kappa_coeff)
    ns = [n1, n2, n3]
    err_vn, err_vs = [], []
    for nv in ns:
        vs, vn, _ = _base.solve_single_step(inputs, {"N": nv}, pole_kappa_coeff=pole_kappa_coeff)
        stride = (N_ref_exact - 1) // (nv - 1)
        err_vn.append(float(np.max(np.abs(vn - vn_ref[::stride]))))
        err_vs.append(float(np.max(np.abs(vs - vs_ref[::stride]))))

    logN = np.log(np.array(ns, dtype=float) - 1.0)
    p_vn = float(-np.polyfit(logN, np.log(np.maximum(err_vn, 1e-16)), 1)[0])
    p_vs = float(-np.polyfit(logN, np.log(np.maximum(err_vs, 1e-16)), 1)[0])
    return {"order_vn": p_vn, "order_vs": p_vs}


def run(inputs, settings):
    return compute_orders(inputs, settings, pole_kappa_coeff=8.0)

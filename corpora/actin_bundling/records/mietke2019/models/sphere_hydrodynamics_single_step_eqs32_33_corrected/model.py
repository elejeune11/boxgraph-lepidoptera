import importlib.util
from pathlib import Path
import numpy as np

_base_path = Path(__file__).resolve().parent.parent / "sphere_hydrodynamics_single_step_eqs32_33" / "model.py"
_spec = importlib.util.spec_from_file_location("_sphere_hydro_base", _base_path)
_base = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_base)


def run(inputs, settings):
    vs, vn, p = _base.solve_single_step(inputs, settings, pole_kappa_coeff=16.0 / 3.0)
    return {
        "vn_max_dimless": float(np.max(np.abs(vn))),
        "vs_max_dimless": float(np.max(np.abs(vs))),
        "vn_pole_dimless": float(vn[0]),
        "p_dimless": float(p),
    }

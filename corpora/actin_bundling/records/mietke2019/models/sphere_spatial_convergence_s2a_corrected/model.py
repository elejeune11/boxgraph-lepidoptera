import importlib.util
from pathlib import Path

_base_path = Path(__file__).resolve().parent.parent / "sphere_spatial_convergence_s2a" / "model.py"
_spec = importlib.util.spec_from_file_location("_sphere_s2a_base", _base_path)
_base = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_base)


def run(inputs, settings):
    return _base.compute_orders(inputs, settings, pole_kappa_coeff=16.0 / 3.0)

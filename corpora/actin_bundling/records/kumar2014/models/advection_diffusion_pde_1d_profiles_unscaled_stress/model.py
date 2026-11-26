import importlib.util
from pathlib import Path

_base_path = Path(__file__).resolve().parent.parent / "advection_diffusion_pde_1d_profiles" / "model.py"
_spec = importlib.util.spec_from_file_location("_ad_pde_1d_profiles", _base_path)
_base = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_base)


def run(inputs, settings):
    return _base.simulate_profiles(inputs, settings, report_unscaled=True)

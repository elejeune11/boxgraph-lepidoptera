import importlib.util
from pathlib import Path

_base_path = Path(__file__).resolve().parent.parent / "sphere_mode_l_stability_eqs69_71" / "model.py"
_spec = importlib.util.spec_from_file_location("_sphere_mode_l_base", _base_path)
_base = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_base)


def run(inputs, settings):
    return _base._solve_mode(inputs, use_corrected_diffusion=True)

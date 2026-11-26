"""Infinite-pool cylinder PDEs with mass-conserving BCs and 1D Cartesian membrane diffusion returning curves."""

import importlib.util
from pathlib import Path

_BASE_PATH = Path(__file__).resolve().parent.parent / "kinetic_pde_cylinder_eq6_16_curves" / "model.py"
_SPEC = importlib.util.spec_from_file_location("_durre_cyl_curves_base", _BASE_PATH)
_BASE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_BASE)


def run(inputs, settings):
    return _BASE._solve_cylinder_curves(inputs, settings, conserved_bc=True, surf_axisym=False)

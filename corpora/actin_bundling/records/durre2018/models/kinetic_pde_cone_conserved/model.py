"""Axisymmetric finite-pool cone kinetic PDEs with mass-conserving BCs (+k_2*V* in Eq. 14, +k_4*A* in Eq. 15)."""

import importlib.util
from pathlib import Path

_BASE_PATH = Path(__file__).resolve().parent.parent / "kinetic_pde_cone_eq6_16" / "model.py"
_SPEC = importlib.util.spec_from_file_location("_durre_cone_base", _BASE_PATH)
_BASE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_BASE)


def run(inputs, settings):
    return _BASE._solve_cone(inputs, settings, conserved_bc=True, surf_axisym=True)

"""Sign-corrected RK4 dimer solver (+delta_{i1}) restoring exact force and momentum conservation."""

from pathlib import Path
import importlib.util

_BASE_PATH = Path(__file__).resolve().parent.parent / "cgmd_dimer_rk4" / "model.py"
_SPEC = importlib.util.spec_from_file_location("_cgmd_dimer_rk4_base", _BASE_PATH)
_MOD = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_MOD)


def run(inputs, settings):
    return _MOD.solve_dimer(inputs, settings, plus_sign=True)

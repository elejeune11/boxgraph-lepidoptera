import importlib.util
from pathlib import Path

_base_path = Path(__file__).resolve().parent.parent / "zero_friction_ode_branch_l2" / "model.py"
_spec = importlib.util.spec_from_file_location("_xi0_base", _base_path)
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)


def run(inputs, settings):
    P = float(inputs["alpha2_phi"])
    i_mode = int(round(inputs["i"]))
    x_targets = [float(v) for v in inputs["target_branch"]["x"]]
    M = int(round(settings.get("M_modes", 40)))
    N = int(round(settings.get("N_quad", 512)))
    return _mod.continue_xi0_ode(P, i_mode, x_targets, "maxrho", M=M, N=N)

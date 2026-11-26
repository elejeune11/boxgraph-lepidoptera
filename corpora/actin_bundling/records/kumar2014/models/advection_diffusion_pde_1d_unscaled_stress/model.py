from pathlib import Path
import importlib.util

_BASE_PATH = Path(__file__).resolve().parent.parent / "advection_diffusion_pde_1d" / "model.py"
_SPEC = importlib.util.spec_from_file_location("_ad1d_base", _BASE_PATH)
_MOD = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_MOD)


def run(inputs, settings):
    res = _MOD.simulate_1d(inputs, settings)
    scale_f = res.pop("scale_f")
    res["f_max"] = res["f_max"] / scale_f
    res["f_min"] = res["f_min"] / scale_f
    return res

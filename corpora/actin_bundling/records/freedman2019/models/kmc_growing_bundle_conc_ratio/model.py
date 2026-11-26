import importlib.util
from pathlib import Path

_BASE_PATH = Path(__file__).resolve().parent.parent / "kmc_growing_bundle_fig3" / "model.py"
_SPEC = importlib.util.spec_from_file_location("_kmc_base_fig3", _BASE_PATH)
_BASE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_BASE)


def run(inputs, settings):
    base_inputs = {
        "mu_grow": inputs["mu_grow"],
        "mu_on_A": float(inputs["conc_ratio"]) * float(inputs["mu_on_B"]),
        "mu_on_B": inputs["mu_on_B"],
        "L_0": inputs["L_0"],
        "l_b": inputs["l_b"],
        "L_p": inputs["L_p"],
        "delta_l_xl": inputs["delta_l_xl"],
        "n_d": inputs["n_d"],
        "l_d": inputs["l_d"],
    }
    return _BASE.run(base_inputs, settings)

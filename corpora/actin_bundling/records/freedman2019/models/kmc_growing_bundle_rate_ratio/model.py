import importlib.util
from pathlib import Path

_BASE_PATH = Path(__file__).resolve().parent.parent / "kmc_growing_bundle_fig3" / "model.py"
_SPEC = importlib.util.spec_from_file_location("_kmc_base_fig3_rr", _BASE_PATH)
_BASE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_BASE)


def run(inputs, settings):
    # Gillespie transition probabilities depend only on r_grow / mu_on = rate_ratio
    l_b = float(inputs["l_b"])
    mu_on = 1.0
    mu_grow = float(inputs["rate_ratio"]) * l_b * mu_on
    base_inputs = {
        "mu_grow": mu_grow,
        "mu_on_A": mu_on,
        "mu_on_B": mu_on,
        "L_0": inputs["L_0"],
        "l_b": l_b,
        "L_p": inputs["L_p"],
        "delta_l_xl": inputs["delta_l_xl"],
        "n_d": inputs["n_d"],
        "l_d": inputs["l_d"],
    }
    return _BASE.run(base_inputs, settings)

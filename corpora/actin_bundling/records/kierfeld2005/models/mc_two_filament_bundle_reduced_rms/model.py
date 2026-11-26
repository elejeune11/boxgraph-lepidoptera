"""Evaluates RMS sqrt(<(|Delta z_12| - l_r)^2>) from two-filament Metropolis MC at reduced potential strength."""
import importlib.util
from pathlib import Path

_base_path = Path(__file__).resolve().parent.parent / "mc_two_filament_bundle" / "model.py"
_spec = importlib.util.spec_from_file_location("_mc_pair_base", _base_path)
_base = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_base)


def run(inputs, settings):
    W_abs = inputs["W_bar_c"] / (1.0 - inputs["reduced_W"])
    _, _, dz_rms = _base.simulate_pair(
        W_abs=W_abs,
        L_p=inputs["L_p"],
        L=inputs["L"],
        l_a=inputs["l_a"],
        l_r=inputs["l_r"],
        n_equil=settings.get("n_equil", 3000),
        n_meas=settings.get("n_meas", 4000),
        seed=settings.get("seed", 1),
    )
    return {"mean_separation": float(dz_rms)}

"""Diagnostic variant of notebook_droplet_ring_slsqp testing the +2 filament-index offset
(Nfil_eff = Nfil + 2) while invoking the exact live notebook functions from VASPDroplet_simulate.ipynb."""

import importlib.util
from pathlib import Path


def _load_base():
    base_path = Path(__file__).resolve().parents[1] / "notebook_droplet_ring_slsqp" / "model.py"
    spec = importlib.util.spec_from_file_location("_graham_base_slsqp", base_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def run(inputs, settings):
    base = _load_base()
    mod_inputs = dict(inputs)
    t_ring = float(inputs["T_ring"])
    rcyl_init = float(inputs["rcyl_init"])
    r_fil = float(inputs["R_fil"])
    n_base = 1.0 if abs(t_ring - rcyl_init) < 1e-6 else (t_ring / r_fil)
    mod_inputs["T_ring"] = (n_base + 2.0) * r_fil
    return base.run(mod_inputs, settings)

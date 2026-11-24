"""BoxGraph software adapter auditing ModelC_2D.grids() and t=0 initialization in ModelC_Cubic.py."""

import contextlib
import importlib.util
import io
import os
import sys
import types
from pathlib import Path

import numpy as np


def _make_h5_shim(store):
    class _MemH5File:
        def __init__(self, path, mode="r"):
            self.path = path
            self.mode = mode

        def create_dataset(self, name, data):
            store[name] = data

        def close(self):
            pass

    mod = types.ModuleType("h5py")
    mod.File = _MemH5File
    return mod


def _load_modelc_module(store):
    shim = _make_h5_shim(store)
    if "h5py" not in sys.modules:
        sys.modules["h5py"] = shim
    src_path = Path(os.environ["BOXGRAPH_SOURCE"]) / "ModelC_Cubic.py"
    spec = importlib.util.spec_from_file_location("_gubbala_ModelC_Cubic_grids", src_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.h5py = shim
    return mod


def run(inputs, settings):
    grid_mult = int(settings.get("grid_mult", 1))
    seed = int(settings.get("seed", 0))

    N = int(round(float(inputs["N"]) * grid_mult))
    L = float(inputs["L"]) * grid_mult

    store = {}
    mc_mod = _load_modelc_module(store)
    mc = mc_mod.ModelC_2D(seed=seed, path="/mem", device="cpu")
    mc.N = N
    mc.L = L
    mc.PHI_AVG = float(inputs["phi0"])
    mc.BETA = float(inputs["beta"])
    mc.GAMMA = float(inputs["gamma"])
    mc.E1 = float(inputs["E1"])
    mc.E3 = float(inputs["E3"])
    mc.CHI = float(inputs["chi"])
    mc.DT = float(inputs["dt_init"])
    mc.TOL = float(inputs["tol"])
    mc.TIME = 0.0
    mc.SKIPS = 1

    with contextlib.redirect_stdout(io.StringIO()):
        mc.run()

    phi0 = np.asarray(store["phi"][0], dtype=float)
    qxx0 = np.asarray(store["qxx"][0], dtype=float)
    qxy0 = np.asarray(store["qxy"][0], dtype=float)

    dx_exported = float(mc.DX)
    dx_fft = float(mc.L / mc.N)
    k_zero_mode = float(mc.K[0, 0])
    dealias_kept_fraction = float(np.mean(mc.ALIASING))
    noise_var_init = float(np.var(qxx0))
    init_qxx_qxy_corr = float(np.corrcoef(qxx0.ravel(), qxy0.ravel())[0, 1])
    init_phi_qxx_corr = float(np.corrcoef((phi0 - mc.PHI_AVG).ravel(), qxx0.ravel())[0, 1])

    return {
        "dx_exported": dx_exported,
        "dx_fft": dx_fft,
        "k_zero_mode": k_zero_mode,
        "dealias_kept_fraction": dealias_kept_fraction,
        "noise_var_init": noise_var_init,
        "init_qxx_qxy_corr": init_qxx_qxy_corr,
        "init_phi_qxx_corr": init_phi_qxx_corr,
    }

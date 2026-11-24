"""BoxGraph software adapter running unmodified ModelC_2D.run() natively on CPU with an in-memory h5py shim."""

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
    spec = importlib.util.spec_from_file_location("_gubbala_ModelC_Cubic_native", src_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.h5py = shim
    return mod


def _build_R_fn(N, L, k_cut_ratio):
    dx = L / N
    kx = 2.0 * np.pi * np.fft.fftfreq(N, d=dx)
    ky = 2.0 * np.pi * np.fft.fftfreq(N, d=dx)
    KX, KY = np.meshgrid(kx, ky)
    K = np.sqrt(KX**2 + KY**2)
    dk = 2.0 * np.pi / L
    nbins = N // 2
    bidx = np.rint(K / dk).astype(int)
    valid = (bidx >= 1) & (bidx < nbins)
    cnt = np.bincount(bidx[valid], minlength=nbins)
    nz = cnt > 0
    k_r = np.arange(nbins) * dk

    def compute_R(phi):
        Pk = np.abs(np.fft.fft2(phi)) ** 2
        Sr = np.bincount(bidx[valid], weights=Pk[valid], minlength=nbins)
        Savg = np.zeros(nbins)
        Savg[nz] = Sr[nz] / cnt[nz]
        kp = int(np.argmax(Savg[1:])) + 1
        sel = (k_r > 0) & (k_r <= k_cut_ratio * k_r[kp])
        den = float(np.sum(k_r[sel] * Savg[sel]))
        return float(np.sum(Savg[sel]) / den) if den > 0.0 else float("nan")

    return compute_R


def run(inputs, settings):
    N = int(settings.get("n_grid", 48))
    L = float(settings.get("L_box", 46.875))
    skips = int(settings.get("skips", 15))
    tol_factor = float(settings.get("tol_factor", 1.0))
    seed = int(settings.get("seed", 0))

    t_fit_start = float(inputs["t_fit_start"])
    t_end = float(inputs["t_end"])
    k_cut_ratio = float(inputs["k_cut_ratio"])

    store = {}
    mc_mod = _load_modelc_module(store)
    mc = mc_mod.ModelC_2D(seed=seed, path="/mem", device="cpu")
    mc.N = N
    mc.L = L
    mc.BETA = float(inputs["beta"])
    mc.GAMMA = float(inputs["gamma"])
    mc.CHI = float(inputs["chi"])
    mc.E1 = float(inputs["E1"])
    mc.E3 = float(inputs["E3"])
    mc.PHI_AVG = float(inputs["phi0"])
    mc.DT = float(inputs["dt_init"])
    mc.TOL = float(inputs["tol"]) * tol_factor
    mc.TIME = t_end
    mc.SKIPS = skips

    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        mc.run()
    log_txt = buf.getvalue()
    blowup_flag = 1.0 if ("Simulation blowup" in log_txt or "picoseconds" in log_txt) else 0.0

    times = np.asarray(store["time"], dtype=float)
    phis = np.asarray(store["phi"], dtype=float)
    qxxs = np.asarray(store["qxx"], dtype=float)
    qxys = np.asarray(store["qxy"], dtype=float)

    if blowup_flag > 0.5 or not np.isfinite(phis[-1]).all():
        return {
            "alpha_fit": float("nan"),
            "R_final": float("nan"),
            "R_growth_ratio": float("nan"),
            "Q_max_final": float("nan"),
            "phi_mean_drift": float("nan"),
            "blowup_flag": 1.0,
        }

    comp_R = _build_R_fn(N, L, k_cut_ratio)
    sel = times >= t_fit_start
    if int(np.sum(sel)) >= 2:
        t_sel = times[sel]
        R_sel = np.array([comp_R(p) for p in phis[sel]], dtype=float)
        alpha_fit = float(np.polyfit(np.log(t_sel), np.log(R_sel), 1)[0])
        R_final = float(R_sel[-1])
        R_growth_ratio = float(R_sel[-1] / R_sel[0]) if R_sel[0] > 0.0 else float("nan")
    else:
        alpha_fit = float("nan")
        R_final = comp_R(phis[-1])
        R_growth_ratio = float("nan")

    Q_max_final = float(np.max(np.hypot(qxxs[-1], qxys[-1])))
    phi_mean_drift = float(abs(np.mean(phis[-1]) - np.mean(phis[0])))

    return {
        "alpha_fit": alpha_fit,
        "R_final": R_final,
        "R_growth_ratio": R_growth_ratio,
        "Q_max_final": Q_max_final,
        "phi_mean_drift": phi_mean_drift,
        "blowup_flag": 0.0,
    }

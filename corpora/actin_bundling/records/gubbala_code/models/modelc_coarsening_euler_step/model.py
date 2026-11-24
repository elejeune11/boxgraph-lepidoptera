"""BoxGraph software adapter driving live ModelC_2D.euler_step() at fixed dt with shared vs independent noise."""

import importlib.util
import os
import sys
import types
from pathlib import Path

import numpy as np


def _load_modelc_module():
    if "h5py" not in sys.modules:
        sys.modules["h5py"] = types.ModuleType("h5py")
    src_path = Path(os.environ["BOXGRAPH_SOURCE"]) / "ModelC_Cubic.py"
    spec = importlib.util.spec_from_file_location("_gubbala_ModelC_Cubic_euler", src_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
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
    N = int(settings.get("n_grid", 64))
    L = float(settings.get("L_box", 62.5))
    dt = float(settings.get("dt", 0.4))
    init_mode = str(settings.get("init_mode", "shared"))
    n_samples = int(settings.get("n_samples", 6))
    seed = int(settings.get("seed", 0))

    phi0_avg = float(inputs["phi0"])
    noise_std = float(inputs["noise_std"])
    t_fit_start = float(inputs["t_fit_start"])
    t_end = float(inputs["t_end"])
    k_cut_ratio = float(inputs["k_cut_ratio"])

    mc_mod = _load_modelc_module()
    mc = mc_mod.ModelC_2D(seed=seed, path="/mem", device="cpu")
    mc.N = N
    mc.L = L
    mc.BETA = float(inputs["beta"])
    mc.GAMMA = float(inputs["gamma"])
    mc.CHI = float(inputs["chi"])
    mc.E1 = float(inputs["E1"])
    mc.E3 = float(inputs["E3"])
    mc.PHI_AVG = phi0_avg
    mc.grids()

    if init_mode == "shared":
        # Exact ModelC_2D.run() lines 148-153 initialization
        noise = mc.prng.normal(loc=0.0, scale=noise_std, size=(N, N))
        phi = noise + phi0_avg
        qxx = noise.copy()
        qxy = noise.copy()
    elif init_mode == "independent":
        # SI p. 17 independent Gaussian perturbations for phi, Qxx, Qxy
        phi = phi0_avg + mc.prng.normal(loc=0.0, scale=noise_std, size=(N, N))
        qxx = mc.prng.normal(loc=0.0, scale=noise_std, size=(N, N))
        qxy = mc.prng.normal(loc=0.0, scale=noise_std, size=(N, N))
    else:
        raise ValueError(f"unknown init_mode: {init_mode!r}")

    phi_mean_init = float(np.mean(phi))
    comp_R = _build_R_fn(N, L, k_cut_ratio)

    nsteps = int(round(t_end / dt))
    sample_times = np.geomspace(t_fit_start, t_end, n_samples)
    sample_steps = {int(round(ts / dt)): float(ts) for ts in sample_times}
    sample_steps[nsteps] = t_end

    ts_rec, Rs_rec = [], []
    for step in range(1, nsteps + 1):
        phi, qxx, qxy = mc.euler_step(phi, qxx, qxy, dt)
        if not np.isfinite(phi).all() or float(np.max(np.abs(phi))) > 1e3:
            return {
                "alpha_fit": float("nan"),
                "R_final": float("nan"),
                "R_growth_ratio": float("nan"),
                "Q_max_final": float("nan"),
                "phi_mean_drift": float("nan"),
                "blowup_flag": 1.0,
            }
        if step in sample_steps:
            ts_rec.append(step * dt)
            Rs_rec.append(comp_R(phi))

    ts_arr = np.asarray(ts_rec, dtype=float)
    Rs_arr = np.asarray(Rs_rec, dtype=float)
    alpha_fit = float(np.polyfit(np.log(ts_arr), np.log(Rs_arr), 1)[0])
    R_final = float(Rs_arr[-1])
    R_growth_ratio = float(Rs_arr[-1] / Rs_arr[0]) if Rs_arr[0] > 0.0 else float("nan")
    Q_max_final = float(np.max(np.hypot(qxx, qxy)))
    phi_mean_drift = float(abs(float(np.mean(phi)) - phi_mean_init))

    return {
        "alpha_fit": alpha_fit,
        "R_final": R_final,
        "R_growth_ratio": R_growth_ratio,
        "Q_max_final": Q_max_final,
        "phi_mean_drift": phi_mean_drift,
        "blowup_flag": 0.0,
    }

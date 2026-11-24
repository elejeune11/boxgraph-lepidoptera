"""BoxGraph software adapter auditing live ModelC_2D.run() adaptive stepping, T += DT_new shift, and mass drift."""

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
    spec = importlib.util.spec_from_file_location("_gubbala_ModelC_Cubic_step", src_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.h5py = shim
    return mod


def _compute_R(phi, N, L, k_cut_ratio):
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
    Pk = np.abs(np.fft.fft2(phi)) ** 2
    Sr = np.bincount(bidx[valid], weights=Pk[valid], minlength=nbins)
    Savg = np.zeros(nbins)
    Savg[nz] = Sr[nz] / cnt[nz]
    kp = int(np.argmax(Savg[1:])) + 1
    sel = (k_r > 0) & (k_r <= k_cut_ratio * k_r[kp])
    den = float(np.sum(k_r[sel] * Savg[sel]))
    return float(np.sum(Savg[sel]) / den) if den > 0.0 else float("nan")


def run(inputs, settings):
    N = int(settings.get("n_grid", 32))
    L = float(settings.get("L_box", 31.25))
    skips = int(settings.get("skips", 1))
    tol_factor = float(settings.get("tol_factor", 1.0))
    seed = int(settings.get("seed", 0))

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
    mc.TIME = float(inputs["t_end"])
    mc.SKIPS = skips

    # Non-intrusive pass-through logger recording the exact dt used on each full euler_step
    dt_used = []
    orig_euler = mc.euler_step
    call_count = [0]

    def logged_euler(phi_, qxx_, qxy_, dt):
        call_count[0] += 1
        if call_count[0] % 3 == 0:
            dt_used.append(float(dt))
        return orig_euler(phi_, qxx_, qxy_, dt)

    mc.euler_step = logged_euler

    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        mc.run()
    log_txt = buf.getvalue()
    blowup_flag = 1.0 if ("Simulation blowup" in log_txt or "picoseconds" in log_txt) else 0.0

    times = np.asarray(store["time"], dtype=float)
    phis = np.asarray(store["phi"], dtype=float)

    # With skips=1, frame k is saved before step k executes; true integrated time at frame k is sum_{j<k} dt_used[j]
    n_frames = len(times)
    true_times = np.concatenate([[0.0], np.cumsum(dt_used)])[:n_frames]

    reported_time_step1 = float(times[1]) if n_frames > 1 else float("nan")
    true_time_step1 = float(true_times[1]) if n_frames > 1 else float("nan")
    max_time_shift = float(np.max(np.abs(times - true_times))) if n_frames > 0 else float("nan")
    phi_mean_drift = float(abs(np.mean(phis[-1]) - np.mean(phis[0])))
    R_final = _compute_R(phis[-1], N, L, float(inputs["k_cut_ratio"]))

    return {
        "reported_time_step1": reported_time_step1,
        "true_time_step1": true_time_step1,
        "max_time_shift": max_time_shift,
        "final_dt": float(mc.DT),
        "phi_mean_drift": phi_mean_drift,
        "R_final": R_final,
        "blowup_flag": blowup_flag,
    }

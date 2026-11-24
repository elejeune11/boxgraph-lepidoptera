"""BoxGraph software adapter probing the linearized operator of live ModelC_2D.euler_step()."""

import importlib.util
import math
import os
import sys
import types
from pathlib import Path

import numpy as np


def _load_modelc_module():
    if "h5py" not in sys.modules:
        sys.modules["h5py"] = types.ModuleType("h5py")
    src_path = Path(os.environ["BOXGRAPH_SOURCE"]) / "ModelC_Cubic.py"
    spec = importlib.util.spec_from_file_location("_gubbala_ModelC_Cubic_lsa", src_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _measure_mode(mc, phi_b, qxx_b, qxy_b, X, Y, kx, ky, dt, eps):
    w = np.cos(kx * X + ky * Y)
    wn = float(np.sum(w * w))

    pp, qxp, _ = mc.euler_step(phi_b + eps * w, qxx_b, qxy_b, dt)
    pm, qxm, _ = mc.euler_step(phi_b - eps * w, qxx_b, qxy_b, dt)
    w_phi = -float(np.sum(((pp - pm) / (2.0 * eps) - w) / dt * w)) / wn
    c21 = float(np.sum(((qxp - qxm) / (2.0 * eps)) / dt * w)) / wn

    pp, qxp, _ = mc.euler_step(phi_b, qxx_b + eps * w, qxy_b, dt)
    pm, qxm, _ = mc.euler_step(phi_b, qxx_b - eps * w, qxy_b, dt)
    c12 = float(np.sum(((pp - pm) / (2.0 * eps)) / dt * w)) / wn
    w_xx = -float(np.sum(((qxp - qxm) / (2.0 * eps) - w) / dt * w)) / wn

    _, _, qyp = mc.euler_step(phi_b, qxx_b, qxy_b + eps * w, dt)
    _, _, qym = mc.euler_step(phi_b, qxx_b, qxy_b - eps * w, dt)
    w_xy = -float(np.sum(((qyp - qym) / (2.0 * eps) - w) / dt * w)) / wn
    return w_phi, c12, c21, w_xx, w_xy


def _extract_k4_coeff(mc, phi_b, qxx_b, qxy_b, X, Y, L, m1, m2, axis, dt, eps):
    """Extract the k^4 coefficient of omega_phi/(gamma*k^2) along longitudinal ('x') or transverse ('y') axis."""
    vals = []
    k2s = []
    for m_idx in (m1, m2):
        k = 2.0 * math.pi * m_idx / L
        kx, ky = (k, 0.0) if axis == "x" else (0.0, k)
        w = np.cos(kx * X + ky * Y)
        wn = float(np.sum(w * w))
        pp, _, _ = mc.euler_step(phi_b + eps * w, qxx_b, qxy_b, dt)
        pm, _, _ = mc.euler_step(phi_b - eps * w, qxx_b, qxy_b, dt)
        w_phi = -float(np.sum(((pp - pm) / (2.0 * eps) - w) / dt * w)) / wn
        vals.append(w_phi / (mc.GAMMA * k * k))
        k2s.append(k * k)
    return (vals[1] - vals[0]) / (k2s[1] - k2s[0])


def _bisect_chi_bound(mc, phi_b, qxx_b, qxy_b, X, Y, L, m1, m2, axis, lo, hi, dt, eps):
    orig_chi = mc.CHI
    try:
        mc.CHI = lo
        f_lo = _extract_k4_coeff(mc, phi_b, qxx_b, qxy_b, X, Y, L, m1, m2, axis, dt, eps)
        mc.CHI = hi
        f_hi = _extract_k4_coeff(mc, phi_b, qxx_b, qxy_b, X, Y, L, m1, m2, axis, dt, eps)
        if f_lo * f_hi > 0.0:
            return float("nan")
        for _ in range(36):
            mid = 0.5 * (lo + hi)
            mc.CHI = mid
            f_mid = _extract_k4_coeff(mc, phi_b, qxx_b, qxy_b, X, Y, L, m1, m2, axis, dt, eps)
            if f_lo * f_mid <= 0.0:
                hi, f_hi = mid, f_mid
            else:
                lo, f_lo = mid, f_mid
        return 0.5 * (lo + hi)
    finally:
        mc.CHI = orig_chi


def run(inputs, settings):
    beta = float(inputs["beta"])
    gamma = float(inputs["gamma"])
    chi = float(inputs["chi"])
    E1 = float(inputs["E1"])
    E3 = float(inputs["E3"])

    N = int(settings.get("n_grid", 32))
    L = float(settings.get("L_box", 32.0))
    m_mode = int(settings.get("mode_m", 2))
    dt = float(settings.get("dt_probe", 1.0e-4))
    eps = float(settings.get("eps_fd", 1.0e-4))

    mc_mod = _load_modelc_module()
    mc = mc_mod.ModelC_2D(seed=0, path="/mem", device="cpu")
    mc.N = N
    mc.L = L
    mc.BETA = beta
    mc.GAMMA = gamma
    mc.CHI = chi
    mc.E1 = E1
    mc.E3 = E3
    mc.grids()

    dx = L / N
    x = np.arange(N) * dx
    X, Y = np.meshgrid(x, x)

    S = 1.0 / math.sqrt(beta)
    phi_b = np.ones((N, N), dtype=float)
    qxx_b = np.full((N, N), S, dtype=float)
    qxy_b = np.zeros((N, N), dtype=float)

    k_val = 2.0 * math.pi * m_mode / L
    k2 = k_val * k_val

    w_phi_l, c12_l, c21_l, w_xx_l, w_xy_l = _measure_mode(mc, phi_b, qxx_b, qxy_b, X, Y, k_val, 0.0, dt, eps)
    w_phi_t, _, _, w_xx_t, _ = _measure_mode(mc, phi_b, qxx_b, qxy_b, X, Y, 0.0, k_val, dt, eps)

    # Under SI Eq. (S28): omega_xx_long - omega_xx_trans = 2 * (E3 / sqrt(beta)) * k^2
    E3_eff_linear = (w_xx_l - w_xx_t) * math.sqrt(beta) / (2.0 * k2)
    # Under Eq. (9): omega_phi_long - omega_phi_trans = 2 * gamma * chi * k^4 (code drops 1/S -> chi/sqrt(beta))
    chi_eff_linear = (w_phi_l - w_phi_t) / (2.0 * gamma * k2 * k2)

    # Live high-k stability bounds on chi from ModelC_2D.euler_step() (k^4 coeff = 0)
    # Use E3=0 temporarily for speed and exactness since chi-stability depends only on phi equation
    mc.E3 = 0.0
    chi_min_code = _bisect_chi_bound(mc, phi_b, qxx_b, qxy_b, X, Y, L, m_mode, 2 * m_mode, "x", -10.0, 0.0, dt, eps)
    chi_max_code = _bisect_chi_bound(mc, phi_b, qxx_b, qxy_b, X, Y, L, m_mode, 2 * m_mode, "y", 0.0, 10.0, dt, eps)
    mc.E3 = E3

    # k = 0 uniform perturbation mass-conservation leak rate from K[0,0] = 1/N^2
    dt_k0 = 1.0e-2
    ones = np.ones((N, N), dtype=float)
    pp0, _, _ = mc.euler_step(phi_b + eps * ones, qxx_b, qxy_b, dt_k0)
    pm0, _, _ = mc.euler_step(phi_b - eps * ones, qxx_b, qxy_b, dt_k0)
    mass_leak_rate = -float(np.mean(((pp0 - pm0) / (2.0 * eps) - ones) / dt_k0))

    return {
        "omega_phi_long": float(w_phi_l),
        "omega_phi_trans": float(w_phi_t),
        "c12_phi_from_qxx": float(c12_l),
        "c21_qxx_from_phi": float(c21_l),
        "omega_xx_long": float(w_xx_l),
        "omega_xx_trans": float(w_xx_t),
        "omega_xy_long": float(w_xy_l),
        "E3_eff_linear": float(E3_eff_linear),
        "chi_eff_linear": float(chi_eff_linear),
        "chi_min_code": float(chi_min_code),
        "chi_max_code": float(chi_max_code),
        "mass_leak_rate": float(mass_leak_rate),
    }

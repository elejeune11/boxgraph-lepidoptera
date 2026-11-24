"""BoxGraph software adapter auditing ModelC_2D.cubic() directly from ModelC_Cubic.py."""

import importlib.util
import math
import os
import sys
import types
from pathlib import Path

import numpy as np


def _load_modelc_module():
    if "h5py" not in sys.modules:
        mod_h5 = types.ModuleType("h5py")
        sys.modules["h5py"] = mod_h5
    src_path = Path(os.environ["BOXGRAPH_SOURCE"]) / "ModelC_Cubic.py"
    spec = importlib.util.spec_from_file_location("_gubbala_ModelC_Cubic_cubic", src_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def run(inputs, settings):
    E3 = float(inputs["E3"])
    beta = float(inputs["beta"])
    N = int(settings.get("n_grid", 32))
    L = float(settings.get("L_box", 32.0))
    eps = float(settings.get("eps_fd", 1.0e-6))

    mc_mod = _load_modelc_module()
    mc = mc_mod.ModelC_2D(seed=0, path="/mem", device="cpu")
    mc.N = N
    mc.L = L
    mc.E3 = E3
    mc.BETA = beta
    mc.grids()

    dx = L / N
    x = np.arange(N) * dx
    X, Y = np.meshgrid(x, x)  # indexing='xy' matching mc.KX, mc.KY

    # Smooth periodic 2D test field (same wave form as gubbala2025/cubic_elastic_variational_sign_eq10)
    Qxx = 0.4 * np.sin(2.0 * np.pi * X / L) * np.cos(2.0 * np.pi * Y / L) + 0.3 * np.cos(4.0 * np.pi * X / L)
    Qxy = 0.3 * np.cos(2.0 * np.pi * X / L) * np.sin(2.0 * np.pi * Y / L) - 0.2 * np.sin(2.0 * np.pi * Y / L)

    qxx_fft = mc.fft2(Qxx)
    qxy_fft = mc.fft2(Qxy)

    # 1. Live ModelC_2D.cubic() call
    fxx_code, fxy_code = mc.cubic(Qxx, qxx_fft, Qxy, qxy_fft)

    # 2. Exact Eq. (10) cubic tensor (using Q_fft in the second-derivative term)
    qxx_x, qxx_y = mc.grad(qxx_fft)
    qxy_x, qxy_y = mc.grad(qxy_fft)
    Qx = np.array([[qxx_x, qxy_x], [qxy_x, -qxx_x]])
    Qy = np.array([[qxx_y, qxy_y], [qxy_y, -qxx_y]])
    grad_Q = np.array([Qx, Qy])
    Q = np.array([[Qxx, Qxy], [Qxy, -Qxx]])
    Q_fft = np.array([[qxx_fft, qxy_fft], [qxy_fft, -qxx_fft]])
    K = np.array([mc.KX, mc.KY])

    F_exact = np.zeros_like(Q)
    F_sec_exact = np.zeros_like(Q)
    F_sec_code = np.zeros_like(Q)
    for i in range(2):
        for j in range(2):
            for k in range(2):
                for l in range(2):
                    t12 = grad_Q[i][k][l] * grad_Q[j][k][l] - 2.0 * grad_Q[l][k][l] * grad_Q[k][i][j]
                    t3_ex = -2.0 * Q[k][l] * mc.ifft2(-K[k] * K[l] * Q_fft[i][j]).real
                    t3_cd = -2.0 * Q[k][l] * mc.ifft2(-K[k] * K[l] * Q[i][j]).real
                    F_exact[i][j] += t12 + t3_ex
                    F_sec_exact[i][j] += t3_ex
                    F_sec_code[i][j] += t3_cd

    tr_ex = F_exact[0][0] + F_exact[1][1]
    fxx_exact = F_exact[0][0] - 0.5 * tr_ex
    fxy_exact = F_exact[0][1]

    diff_tot = np.hypot(fxx_code - fxx_exact, fxy_code - fxy_exact)
    norm_tot = np.hypot(fxx_exact, fxy_exact)
    rel_err_vs_eq10 = float(np.linalg.norm(diff_tot) / np.linalg.norm(norm_tot))

    diff_sec = np.hypot(F_sec_code[0][0] - F_sec_exact[0][0], F_sec_code[0][1] - F_sec_exact[0][1])
    norm_sec = np.hypot(F_sec_exact[0][0], F_sec_exact[0][1])
    rel_err_second_deriv = float(np.linalg.norm(diff_sec) / np.linalg.norm(norm_sec))

    # 3. Free-energy dissipation alignment ratio for live euler_step's +0.5 * E3 * cubic()
    def eval_F_E3(qxx, qxy):
        qxk, qyk = mc.fft2(qxx), mc.fft2(qxy)
        dxx_x, dxx_y = mc.grad(qxk)
        dxy_x, dxy_y = mc.grad(qyk)
        gxx = 2.0 * (dxx_x**2 + dxy_x**2)
        gyy = 2.0 * (dxx_y**2 + dxy_y**2)
        gxy = 2.0 * (dxx_x * dxx_y + dxy_x * dxy_y)
        return float(np.sum(0.5 * E3 * (qxx * (gxx - gyy) + 2.0 * qxy * gxy)) * dx * dx)

    rxx_code = 0.5 * E3 * fxx_code
    rxy_code = 0.5 * E3 * fxy_code
    F0 = eval_F_E3(Qxx, Qxy)
    dF_dt = (eval_F_E3(Qxx + eps * rxx_code, Qxy + eps * rxy_code) - F0) / eps
    expected = -2.0 * float(np.sum(rxx_code**2 + rxy_code**2) * dx * dx)
    dissipation_ratio_code = float(dF_dt / expected) if expected != 0.0 else float("nan")

    # 4. Uniform constant nematic field residual (exact mathematical value = 0 everywhere)
    S = 1.0 / math.sqrt(beta)
    qxx_u = np.full((N, N), S, dtype=float)
    qxy_u = np.zeros((N, N), dtype=float)
    fxx_u, fxy_u = mc.cubic(qxx_u, mc.fft2(qxx_u), qxy_u, mc.fft2(qxy_u))
    uniform_residual_max = float(max(np.max(np.abs(0.5 * E3 * fxx_u)), np.max(np.abs(0.5 * E3 * fxy_u))))

    return {
        "dissipation_ratio_code": dissipation_ratio_code,
        "rel_err_vs_eq10": rel_err_vs_eq10,
        "rel_err_second_deriv": rel_err_second_deriv,
        "uniform_residual_max": uniform_residual_max,
    }

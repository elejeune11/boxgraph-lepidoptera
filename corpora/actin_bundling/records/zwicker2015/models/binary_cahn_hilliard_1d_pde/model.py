import math
import numpy as np
from scipy.optimize import brentq


def _integrate_1d(kf_h, kb_h, L, nx, dt, T, phi_init):
    k = 2.0 * np.pi * np.fft.rfftfreq(nx, d=L / nx)
    k2 = k**2
    k4 = k**4
    # Stabilized linear operator: L_op = -2*k2 - 0.25*k4 - (kf_h + kb_h)
    # Nonlinear residual: d_xx(phi - 3*phi^2 + 2*phi^3) - (-2*d_xx phi) + kf_h = d_xx(-phi - 3*phi^2 + 2*phi^3) + kf_h
    S = 2.0
    L_op = -S * k2 - 0.25 * k4 - (kf_h + kb_h)
    # Second-order IMEX (Crank-Nicolson linear + midpoint/Heun nonlinear)
    lhs = 1.0 - 0.5 * dt * L_op
    rhs_lin = 1.0 + 0.5 * dt * L_op
    nsteps = max(1, int(round(T / dt)))
    h = T / nsteps

    phi = phi_init.copy()
    phi_hat = np.fft.rfft(phi)
    phi_min = float(np.min(phi))

    def nl_hat_fn(p):
        nl = (1.0 - S) * p - 3.0 * p**2 + 2.0 * p**3
        nh = -k2 * np.fft.rfft(nl)
        nh[0] += kf_h * nx
        return nh

    # Recompute operators for exact h
    lhs = 1.0 - 0.5 * h * L_op
    rhs_lin = 1.0 + 0.5 * h * L_op
    lhs_e = 1.0 - h * L_op

    for _ in range(nsteps):
        nh0 = nl_hat_fn(phi)
        # Predictor (semi-implicit Euler)
        phi_pred_hat = (phi_hat + h * nh0) / lhs_e
        phi_pred = np.fft.irfft(phi_pred_hat, n=nx)
        nh1 = nl_hat_fn(phi_pred)
        # Corrector (Crank-Nicolson + trapezoidal nonlinear)
        phi_hat = (rhs_lin * phi_hat + 0.5 * h * (nh0 + nh1)) / lhs
        phi = np.fft.irfft(phi_hat, n=nx)
        pmin = float(np.min(phi))
        if pmin < phi_min:
            phi_min = pmin
    return phi, phi_min


def _measure_mode_rate(kf_h, kb_h, qw, L, nx, dt, eps, T):
    # Pick integer mode closest to qw * L / (2*pi), and adjust L_eff so q_mode == qw exactly
    m = max(1, int(round(qw * L / (2.0 * math.pi))))
    L_eff = 2.0 * math.pi * m / qw
    x = np.linspace(0.0, L_eff, nx, endpoint=False)
    phi0 = kf_h / (kf_h + kb_h)
    phi_init = phi0 + eps * np.cos(qw * x)
    phi_end, phi_min = _integrate_1d(kf_h, kb_h, L_eff, nx, dt, T, phi_init)
    amp_end = 2.0 * abs(np.fft.rfft(phi_end - phi0)[m]) / nx
    lam = math.log(max(amp_end, 1e-300) / eps) / T
    return lam, phi_min


def run(inputs, settings):
    kf_h = inputs["kf_norm"]
    kb_h = inputs["kb_norm"]
    qw = inputs["qw"]
    L = inputs["L_over_w"]
    eps = inputs["eps"]
    T = inputs["T_k0"]
    nx = int(round(settings.get("nx", 128)))
    dt = float(settings.get("dt", 0.02))

    # 1. Measure linear mode growth rate at requested (kf_norm, kb_norm, qw)
    lam_mode, phi_min_mode = _measure_mode_rate(kf_h, kb_h, qw, L, nx, dt, eps, T)

    # 2. Directly locate the symmetric (kf = kb) homogeneous stability cutoff at peak wavenumber qw = 1
    def lam_sym_at_k(k_sym):
        val, _ = _measure_mode_rate(k_sym, k_sym, 1.0, L, nx, dt, eps, min(T, 2.0))
        return val

    k_crit = brentq(lam_sym_at_k, 0.02, 0.6, xtol=1e-6, rtol=1e-6)
    max_rate_sum_norm = 2.0 * k_crit  # (kf + kb) / k0 at kf = kb

    # 3. Passive flat-interface equilibrium relaxation to verify w and gamma/(b*w) = 1/12
    L_if = max(L, 30.0)
    x_if = np.linspace(0.0, L_if, nx, endpoint=False)
    phi_init_if = 0.5 * (np.tanh((x_if - 0.25 * L_if) / 1.4) - np.tanh((x_if - 0.75 * L_if) / 1.4))
    phi_if, phi_min_if = _integrate_1d(0.0, 0.0, L_if, nx, dt, 25.0, phi_init_if)
    k_if = 2.0 * np.pi * np.fft.rfftfreq(nx, d=L_if / nx)
    dphi_dx = np.fft.irfft(1j * k_if * np.fft.rfft(phi_if), n=nx)
    f_dens = 0.5 * (phi_if**2) * ((1.0 - phi_if) ** 2) + 0.125 * (dphi_dx**2)
    gamma_over_bw = float(np.sum(f_dens) * (L_if / nx) / 2.0)
    w_num_over_w = float(1.0 / (2.0 * np.max(np.abs(dphi_dx))))

    return {
        "lam_mode_norm": float(lam_mode),
        "phi_min_mode": float(phi_min_mode),
        "max_rate_sum_norm": float(max_rate_sum_norm),
        "w_num_over_w": w_num_over_w,
        "gamma_over_bw": gamma_over_bw,
        "phi_min_interface": float(phi_min_if),
    }

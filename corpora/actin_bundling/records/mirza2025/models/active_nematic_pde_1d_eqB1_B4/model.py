import math
import numpy as np
from scipy.linalg import solve_banded
from scipy.optimize import minimize_scalar
from scipy.signal import find_peaks

def run(inputs, settings):
    kd = float(inputs["k_d"])
    Lf = float(inputs["L"])
    c0 = float(inputs["c_0"])
    b = float(inputs["b"])
    er = float(inputs["eta_rot"])
    beta = float(inputs["beta"])
    lo = float(inputs["lambda_o"])
    kap = float(inputs["kappa"])
    ratio = float(inputs["activity_ratio"])
    Ldom = float(inputs["L_domain"])

    two_a = c0 + lo
    G = c0 if c0 >= 0.0 else -2.0 * c0
    q0 = 0.0 if c0 >= 0.0 else -0.5 * math.sqrt(-c0 / b)

    def lam_fn(x):
        if x <= 0.0:
            return 1e30
        num = (1.0 / x + 4.0) * (x + kd) * (2.0 * G + 2.0 * Lf * x)
        den = (1.0 + kap * q0) * (2.0 * G + 2.0 * Lf * x) + kap * (beta * (x + kd) + 2.0 * lo * q0)
        return num / den if den > 0.0 else 1e30

    xs = np.geomspace(1e-4, 1e3, 2000)
    vals = np.array([lam_fn(float(x)) for x in xs])
    i0 = int(np.argmin(vals))
    res = minimize_scalar(
        lam_fn,
        bounds=(float(xs[max(0, i0 - 3)]), float(xs[min(len(xs) - 1, i0 + 3)])),
        method="bounded",
        options={"xatol": 1e-12},
    )
    lam_c = float(res.fun)
    nu_c = math.sqrt(float(res.x))
    ell_pat = 2.0 * math.pi / nu_c
    lam = ratio * lam_c

    N = int(round(float(settings.get("N_grid", 64))))
    dt = float(settings.get("dt", 0.02))
    t_max = float(settings.get("t_max", 120.0))
    init_mode = int(settings.get("init_mode", 4))
    seed = int(settings.get("seed", 42))

    dx = Ldom / N
    x = np.linspace(0.0, Ldom, N, endpoint=False)
    if init_mode > 0:
        rho = 1.0 + 0.05 * np.cos(2.0 * np.pi * init_mode * x / Ldom)
        q = q0 - 0.05 * np.cos(2.0 * np.pi * init_mode * x / Ldom)
    else:
        rng = np.random.default_rng(seed)
        rho = np.ones(N)
        q = np.full(N, q0)
        # Smooth low-mode deterministic perturbation invariant to grid size N
        m_max = max(6, int(round(1.5 * Ldom / ell_pat)))
        for m in range(1, m_max + 1):
            ph1, ph2 = rng.uniform(0.0, 2.0 * np.pi, 2)
            arg = 2.0 * np.pi * m * x / Ldom
            rho += (0.01 / m_max) * np.cos(arg + ph1)
            q += (0.005 / m_max) * np.cos(arg + ph2)

    k = 2.0 * np.pi * np.fft.rfftfreq(N, d=dx)
    k2 = k * k
    mu_eff = 4.0 - (beta * beta) / (2.0 * er)

    def d1(f):
        return (np.roll(f, -1) - np.roll(f, 1)) / (2.0 * dx)

    def d2(f):
        return (np.roll(f, -1) - 2.0 * f + np.roll(f, 1)) / (dx * dx)

    def solve_v(rho_arr, rhs):
        a_face = mu_eff * 0.5 * (rho_arr + np.roll(rho_arr, -1))
        a_m = np.roll(a_face, 1)
        main = rho_arr + (a_face + a_m) / (dx * dx)
        upper = -a_face / (dx * dx)
        lower = -a_m / (dx * dx)
        g = -main[0]
        main_t = main.copy()
        main_t[0] -= g
        main_t[-1] -= (upper[-1] * lower[0]) / g
        ab = np.zeros((3, N))
        ab[0, 1:] = upper[:-1]
        ab[1, :] = main_t
        ab[2, :-1] = lower[1:]
        u = np.zeros(N)
        u[0] = g
        u[-1] = upper[-1]
        v_vec = np.zeros(N)
        v_vec[0] = 1.0
        v_vec[-1] = lower[0] / g
        Z = solve_banded((1, 1), ab, np.column_stack([rhs, u]), overwrite_ab=True, overwrite_b=True)
        y, z = Z[:, 0], Z[:, 1]
        return y - (np.dot(v_vec, y) / (1.0 + np.dot(v_vec, z))) * z

    denom_rho = 1.0 + dt * (k2 + kd)
    denom_q = 1.0 + dt * (Lf / er) * k2
    n_steps = int(round(t_max / dt))

    for _ in range(n_steps):
        rho_x = d1(rho)
        q_x = d1(q)
        q_xx = d2(q)
        H_q = (two_a - rho * lo + 4.0 * b * q * q) * q - Lf * q_xx - (Lf / rho) * q_x * rho_x
        sigma_rhs = rho * (lam * (1.0 + kap * q) - 2.0 * Lf * q_x * q_x - (beta / er) * H_q)
        v = solve_v(rho, d1(sigma_rhs))
        v_x = d1(v)
        nl_rho = -d1(rho * v) + kd
        rho = np.fft.irfft(np.fft.rfft(rho + dt * nl_rho) / denom_rho, n=N)
        rho = np.maximum(rho, 1e-3)
        rhs_q = -v * q_x - (beta / (2.0 * er)) * v_x + (1.0 / er) * (rho * lo * q + (Lf / rho) * q_x * rho_x)
        loc_denom = 1.0 + dt * (two_a + 4.0 * b * q * q) / er
        q_star = (q + dt * rhs_q) / loc_denom
        q = np.fft.irfft(np.fft.rfft(q_star) / denom_q, n=N)

    rho_x = d1(rho)
    q_x = d1(q)
    q_xx = d2(q)
    H_q = (two_a - rho * lo + 4.0 * b * q * q) * q - Lf * q_xx - (Lf / rho) * q_x * rho_x
    sigma_rhs = rho * (lam * (1.0 + kap * q) - 2.0 * Lf * q_x * q_x - (beta / er) * H_q)
    v = solve_v(rho, d1(sigma_rhs))
    v_x = d1(v)
    q_hat = -(beta / (2.0 * er)) * v_x - (1.0 / er) * H_q

    # For kap < 0 (q < 0 on dense bundles), bundle axis || is y (sigma_22) and perp is x (sigma_11)
    sig_perp_act = rho * lam * (1.0 + kap * q)
    sig_par_act = rho * lam * (1.0 - kap * q)
    sig_perp_vis = rho * (4.0 * v_x + beta * q_hat)
    sig_par_vis = rho * (2.0 * v_x - beta * q_hat)
    sig_perp_nem = -2.0 * rho * Lf * q_x * q_x
    sig_perp = sig_perp_act + sig_perp_vis + sig_perp_nem
    sig_par = sig_par_act + sig_par_vis

    i_pk = int(np.argmax(rho))
    # Periodic peak counting via circular padding
    pad = np.concatenate([rho[-2:], rho, rho[:2]])
    pks, _ = find_peaks(pad, height=1.15)
    pks = pks[(pks >= 2) & (pks < N + 2)]
    n_pks = max(1, len(pks))
    w_ratio = (Ldom / n_pks) / ell_pat

    return {
        "rho_max": float(rho[i_pk]),
        "rho_min": float(np.min(rho)),
        "S_max": float(2.0 * np.max(np.abs(q))),
        "sigma_par_peak_over_lam": float(sig_par[i_pk] / lam),
        "sigma_perp_peak_over_lam": float(sig_perp[i_pk] / lam),
        "sigma_diff_peak_over_lam": float((sig_par[i_pk] - sig_perp[i_pk]) / lam),
        "sigma_par_peak_over_lamc": float(sig_par[i_pk] / lam_c),
        "sigma_perp_peak_over_lamc": float(sig_perp[i_pk] / lam_c),
        "sigma_act_par_peak_over_lam": float(sig_par_act[i_pk] / lam),
        "sigma_act_perp_peak_over_lam": float(sig_perp_act[i_pk] / lam),
        "sigma_vis_par_peak_over_lam": float(sig_par_vis[i_pk] / lam),
        "sigma_vis_perp_peak_over_lam": float(sig_perp_vis[i_pk] / lam),
        "wavelength_ratio": float(w_ratio),
    }

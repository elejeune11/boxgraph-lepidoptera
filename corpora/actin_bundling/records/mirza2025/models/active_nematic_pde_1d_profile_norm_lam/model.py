import math
import numpy as np
from scipy.linalg import solve_banded
from scipy.optimize import minimize_scalar


def _solve_1d_profiles(inputs, settings, norm_by_lam: bool):
    kd = float(inputs["k_d"])
    Lf = float(inputs["L"])
    c0 = float(inputs["c_0"])
    b = float(inputs["b"])
    er = float(inputs["eta_rot"])
    beta = float(inputs["beta"])
    lo = float(inputs["lambda_o"])
    kap = float(inputs["kappa"])
    ratio = float(inputs["activity_ratio"])
    Lper = float(inputs["L_period"])
    rho_ref = inputs["rho_ref"]
    xs_eval = np.asarray(rho_ref["x"], dtype=float)
    ys_ref = np.asarray(rho_ref["y"], dtype=float)

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
    lam = ratio * lam_c

    N = int(round(float(settings.get("N_grid", 64))))
    dt = float(settings.get("dt", 0.04))
    t_max = float(settings.get("t_max", 120.0))
    steady_tol = float(settings.get("steady_tol", 1e-7))

    dx = Lper / N
    x_grid = np.linspace(0.0, Lper, N, endpoint=False)
    rho = 1.0 + 0.08 * np.cos(2.0 * np.pi * x_grid / Lper)
    q = q0 - 0.08 * np.cos(2.0 * np.pi * x_grid / Lper)

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
    min_steps = int(round(25.0 / dt))

    for step in range(1, n_steps + 1):
        rho_prev = rho
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
        if step >= min_steps and step % 25 == 0:
            if float(np.max(np.abs(rho - rho_prev))) < steady_tol * dt:
                break

    rho_x = d1(rho)
    q_x = d1(q)
    q_xx = d2(q)
    H_q = (two_a - rho * lo + 4.0 * b * q * q) * q - Lf * q_xx - (Lf / rho) * q_x * rho_x
    sigma_rhs = rho * (lam * (1.0 + kap * q) - 2.0 * Lf * q_x * q_x - (beta / er) * H_q)
    v = solve_v(rho, d1(sigma_rhs))
    v_x = d1(v)
    q_hat = -(beta / (2.0 * er)) * v_x - (1.0 / er) * H_q

    scale = lam if norm_by_lam else lam_c
    sig_perp_act = rho * lam * (1.0 + kap * q) / scale
    sig_par_act = rho * lam * (1.0 - kap * q) / scale
    sig_perp_vis = rho * (4.0 * v_x + beta * q_hat) / scale
    sig_par_vis = rho * (2.0 * v_x - beta * q_hat) / scale
    sig_perp_nem = -2.0 * rho * Lf * q_x * q_x / scale
    sig_perp = sig_perp_act + sig_perp_vis + sig_perp_nem
    sig_par = sig_par_act + sig_par_vis

    # Periodic trigonometric interpolation helper
    def interp_periodic(field_arr, x_pts):
        fh = np.fft.rfft(field_arr) / N
        m = np.arange(len(fh))
        theta = 2.0 * np.pi * (x_pts[:, None] % Lper) / Lper
        weights = np.ones(len(fh))
        weights[1:-1] = 2.0
        if N % 2 == 1:
            weights[-1] = 2.0
        return np.real(np.sum(fh[None, :] * weights[None, :] * np.exp(1j * theta * m[None, :]), axis=1))

    # Align the single arbitrary periodic translation phase x_shift in [0, Lper) to rho_ref
    shifts = np.linspace(0.0, Lper, 280, endpoint=False)
    best_s = 0.0
    best_err = 1e30
    for s in shifts:
        pred = interp_periodic(rho, xs_eval - s)
        err = float(np.sum((pred - ys_ref) ** 2))
        if err < best_err:
            best_err = err
            best_s = float(s)
    # Fine refinement around best_s
    res_s = minimize_scalar(
        lambda s: float(np.sum((interp_periodic(rho, xs_eval - s) - ys_ref) ** 2)),
        bounds=(best_s - 0.15, best_s + 0.15),
        method="bounded",
        options={"xatol": 1e-6},
    )
    x_shift = float(res_s.x)
    x_shifted = xs_eval - x_shift
    x_out = [float(xv) for xv in xs_eval]
    vis_ref = inputs.get("vis_ref", rho_ref)
    xs_vis = np.asarray(vis_ref["x"], dtype=float)
    x_vis_shifted = xs_vis - x_shift
    x_vis_out = [float(xv) for xv in xs_vis]

    return {
        "rho": {"x": x_out, "y": [float(v) for v in interp_periodic(rho, x_shifted)]},
        "sigma_par": {"x": x_out, "y": [float(v) for v in interp_periodic(sig_par, x_shifted)]},
        "sigma_perp": {"x": x_out, "y": [float(v) for v in interp_periodic(sig_perp, x_shifted)]},
        "sigma_act_par": {"x": x_out, "y": [float(v) for v in interp_periodic(sig_par_act, x_shifted)]},
        "sigma_act_perp": {"x": x_out, "y": [float(v) for v in interp_periodic(sig_perp_act, x_shifted)]},
        "sigma_vis_par": {"x": x_vis_out, "y": [float(v) for v in interp_periodic(sig_par_vis, x_vis_shifted)]},
        "sigma_vis_perp": {"x": x_vis_out, "y": [float(v) for v in interp_periodic(sig_perp_vis, x_vis_shifted)]},
    }


def run(inputs, settings):
    return _solve_1d_profiles(inputs, settings, norm_by_lam=True)

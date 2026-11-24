import math
import numpy as np
from scipy.optimize import minimize


def _binom_half_sq(n):
    val = 1.0
    for k in range(n):
        val *= (0.5 - k) / (k + 1)
    return val * val


def _ellipse_perim_and_da(a, R, c_perim, n_arr):
    b = (R * R) / a
    db_da = -(R * R) / (a * a)
    apb = a + b
    amb = a - b
    h = (amb / apb) ** 2
    dh_da = 2.0 * (amb / apb) * (2.0 * b - 2.0 * a * db_da) / (apb * apb)
    hp = h ** n_arr
    S = np.sum(c_perim * hp)
    dS_dh = np.sum(c_perim[1:] * n_arr[1:] * (h ** (n_arr[1:] - 1))) if len(n_arr) > 1 else 0.0
    P = math.pi * apb * S
    dP_da = math.pi * (1.0 + db_da) * S + math.pi * apb * dS_dh * dh_da
    return P, dP_da


def run(inputs, settings):
    # Literally implements printed Eq. S2 (no rest length L_cyl in stretch, +cos in bend)
    # and printed Eqs. S4-S5 (Delta = 0.5*(gamma_AVR - gamma_AVD), gamma(e) = Delta*tanh(f(e)) + Delta + gamma_AVD)
    D_nm = inputs["D_droplet"] * 1000.0
    R = 0.5 * D_nm
    T_ring = inputs["T_ring"]
    delta = inputs["delta"]
    Lp_nm = inputs["L_p"] * 1000.0
    kBT = inputs["kBT"]
    L_cyl = inputs["L_cyl"]
    k_stretch_fil = inputs["k_stretch_fil"]
    gamma_AVR = inputs["gamma_AVR"]
    gamma_AVD = inputs["gamma_AVD"]
    k_e = inputs["k_e"]
    gamma_I = inputs["gamma_I"]
    L_ratio = inputs["L_ratio"]
    noise_amp = inputs["noise_amp"]

    seed = int(settings.get("seed", 42))
    ftol = float(settings.get("ftol", 1e-6))
    n_wet = max(1, int(round(settings.get("n_wet_subseg", 10))))
    n_terms = max(1, int(round(settings.get("n_perim_terms", 7))))
    r_init_frac = float(settings.get("r_init_frac", 0.96))
    maxiter = max(10, int(round(settings.get("maxiter", 300))))

    c_perim = np.array([_binom_half_sq(n) for n in range(n_terms)])
    n_arr = np.arange(n_terms)

    N_fil = T_ring / delta
    k_stretch = N_fil * k_stretch_fil
    k_bend = N_fil * (Lp_nm * kBT / L_cyl)

    # Printed Eq. S4-S5:
    Delta = 0.5 * (gamma_AVR - gamma_AVD)

    N_seg = max(3, int(round(L_ratio * D_nm / L_cyl)))
    N_beads = N_seg + 1
    r_init = r_init_frac * R
    theta = np.linspace(0.0, (N_seg * L_cyl) / r_init, N_beads)
    rng = np.random.default_rng(seed)
    coords0 = np.column_stack([r_init * np.cos(theta), r_init * np.sin(theta)])
    coords0 = coords0 + rng.uniform(0.0, noise_amp, coords0.shape)

    s_pts = np.linspace(0.0, 1.0, n_wet + 1)
    w1 = (1.0 - s_pts)[None, :]
    w2 = s_pts[None, :]
    eps2 = 1e-6  # smooth regularization for zero segment length under printed Eq. S2

    def energy_and_grad(x):
        a = x[0]
        b = (R * R) / a
        db_da = -(R * R) / (a * a)
        coords = x[1:].reshape(N_beads, 2)
        grad_coords = np.zeros_like(coords)
        grad_a = 0.0

        # Printed Eq. S2: r_{b,i} - r_{b,i+1}
        dr = coords[:-1] - coords[1:]
        seg_sq = np.sum(dr * dr, axis=1)
        seg_len = np.sqrt(seg_sq + eps2)
        u = dr / seg_len[:, None]

        # Printed Eq. S2 stretching term: sum k_stretch * ||r_i - r_{i+1}||^2 (no L_cyl)
        E_stretch = k_stretch * np.sum(seg_sq)
        g_dr_s = (2.0 * k_stretch) * dr
        grad_coords[:-1] += g_dr_s
        grad_coords[1:] -= g_dr_s

        # Printed Eq. S2 bending term: + sum k_bend * (u_i . u_{i+1})
        u1 = u[:-1]
        u2 = u[1:]
        dots = np.sum(u1 * u2, axis=1)
        E_bend = k_bend * np.sum(dots)
        g_dr1 = k_bend * (u2 - dots[:, None] * u1) / seg_len[:-1, None]
        g_dr2 = k_bend * (u1 - dots[:, None] * u2) / seg_len[1:, None]
        grad_coords[:-2] += g_dr1
        grad_coords[1:-1] += g_dr2 - g_dr1
        grad_coords[2:] -= g_dr2

        # Printed Eq. S3-S7 wetting term (note dr_fwd = r_{i+1} - r_i = -dr)
        dr_fwd = coords[1:] - coords[:-1]
        pts = coords[:-1, None, :] + s_pts[None, :, None] * dr_fwd[:, None, :]
        e_pts = (pts[:, :, 0] / a) ** 2 + (pts[:, :, 1] / b) ** 2
        e_avg = 0.5 * (e_pts[:, :-1] + e_pts[:, 1:])
        arg = np.clip(k_e * (e_avg - 1.0), -50.0, 50.0)
        th = np.tanh(arg)
        gamma_e = Delta * th + Delta + gamma_AVD
        dL = seg_len[:, None] / n_wet
        E_wet = T_ring * np.sum(gamma_e * dL)

        dEw_dseg = T_ring * np.sum(gamma_e, axis=1) / n_wet
        g_dr_w = (dEw_dseg / seg_len)[:, None] * dr_fwd
        grad_coords[1:] += g_dr_w
        grad_coords[:-1] -= g_dr_w

        sech2 = 1.0 - th * th
        dEw_deavg = T_ring * dL * Delta * k_e * sech2
        dEw_depts = np.zeros_like(e_pts)
        dEw_depts[:, :-1] += 0.5 * dEw_deavg
        dEw_depts[:, 1:] += 0.5 * dEw_deavg

        de_da = -2.0 * (pts[:, :, 0] ** 2) / (a ** 3) - 2.0 * (pts[:, :, 1] ** 2) / (b ** 3) * db_da
        grad_a += np.sum(dEw_depts * de_da)

        de_dx = 2.0 * pts[:, :, 0] / (a * a)
        de_dy = 2.0 * pts[:, :, 1] / (b * b)
        dEw_dpx = dEw_depts * de_dx
        dEw_dpy = dEw_depts * de_dy
        grad_coords[:-1, 0] += np.sum(dEw_dpx * w1, axis=1)
        grad_coords[:-1, 1] += np.sum(dEw_dpy * w1, axis=1)
        grad_coords[1:, 0] += np.sum(dEw_dpx * w2, axis=1)
        grad_coords[1:, 1] += np.sum(dEw_dpy * w2, axis=1)

        P, dP_da = _ellipse_perim_and_da(a, R, c_perim, n_arr)
        E_surf = gamma_I * P
        grad_a += gamma_I * dP_da

        return float(E_stretch + E_bend + E_wet + E_surf), np.concatenate([[grad_a], grad_coords.ravel()])

    x0 = np.concatenate([[R], coords0.ravel()])
    bounds = [(R * 0.25, R * 4.0)] + [(-3.0 * D_nm, 3.0 * D_nm)] * (2 * N_beads)
    res = minimize(
        energy_and_grad,
        x0,
        jac=True,
        method="SLSQP",
        bounds=bounds,
        options={"ftol": ftol, "maxiter": maxiter},
    )
    a_opt = float(res.x[0])
    b_opt = float((R * R) / a_opt)
    ar = max(a_opt / b_opt, b_opt / a_opt)
    return {
        "aspect_ratio": ar,
        "total_energy": float(res.fun),
        "a": max(a_opt, b_opt) / 1000.0,
        "b": min(a_opt, b_opt) / 1000.0,
    }

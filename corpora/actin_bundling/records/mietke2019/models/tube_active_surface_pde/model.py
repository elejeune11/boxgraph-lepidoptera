import math
import numpy as np
from scipy.integrate import solve_ivp


def spectral_periodic_ops(N, L0):
    du = L0 / N
    freqs = np.fft.fftfreq(N, d=du) * 2.0 * math.pi
    idx = np.arange(N)
    F = np.exp(-2j * math.pi * np.outer(idx, idx) / N)
    Finv = np.exp(2j * math.pi * np.outer(idx, idx) / N) / N
    D1 = np.real(Finv @ ((1j * freqs)[:, None] * F))
    D2 = np.real(Finv @ ((-(freqs * freqs))[:, None] * F))
    return du, D1, D2


def solve_hydro_tube(h, r, Css, Cphi, Gamma_u, c, xi, kappa, eta_ratio, m, D1, D2, du):
    N = len(h)
    ep = eta_ratio + 1.0
    em = 1.0 - eta_ratio
    inv_h = 1.0 / h
    Gamma_s = Gamma_u * inv_h
    sin_psi = r * Cphi
    cos_psi = r * Gamma_s
    Ckk = Css + Cphi

    Ds = inv_h[:, None] * D1
    dh_du = D1 @ h
    dCss_ds = Ds @ Css
    dCphi_ds = Gamma_s * (Css - Cphi)
    dGamma_ds = -Cphi * Css - Gamma_s * Gamma_s

    if xi != 0.0:
        cm = c ** m
        fc = cm / (1.0 + cm)
        dfc_ds = Ds @ fc
    else:
        fc = np.zeros(N)
        dfc_ds = np.zeros(N)

    diag_idx = np.arange(N)
    A_n = ep * (Css * Css + Cphi * Cphi) + 2.0 * em * Css * Cphi
    inv_An = 1.0 / A_n
    B_n = (ep * Css + em * Cphi)[:, None] * Ds
    B_n[diag_idx, diag_idx] += (ep * Cphi + em * Css) * Gamma_s

    rhs_n = -(xi * fc + 1.0) * Ckk
    if kappa != 0.0:
        dCkk_du = D1 @ Ckk
        d2Ckk = (inv_h * inv_h) * (D2 @ Ckk) - (dh_du * (inv_h ** 3)) * dCkk_du
        rhs_n = rhs_n + kappa * Ckk * ((Css - Cphi) ** 2) + 2.0 * kappa * (d2Ckk + Gamma_s * inv_h * dCkk_du)

    Vn_vs = -inv_An[:, None] * B_n
    Vn_p = inv_An
    Vn_lz = inv_An * cos_psi
    Vn_rhs = inv_An * rhs_n

    Ds2 = (inv_h * inv_h)[:, None] * D2 - (dh_du * (inv_h ** 3))[:, None] * D1
    M_t_vs = ep * Ds2 + (em + 2.0 * eta_ratio) * Gamma_s[:, None] * Ds
    M_t_vs[diag_idx, diag_idx] += em * dGamma_ds - 2.0 * eta_ratio * (Gamma_s ** 2)

    M_t_vn = (ep * Css + em * Cphi)[:, None] * Ds
    M_t_vn[diag_idx, diag_idx] += ep * dCss_ds + em * dCphi_ds + 2.0 * eta_ratio * Gamma_s * (Css - Cphi)

    M = np.empty((N + 2, N + 2))
    rhs = np.empty(N + 2)
    M[0:N, 0:N] = M_t_vs + M_t_vn @ Vn_vs
    M[0:N, N] = M_t_vn @ Vn_p
    M[0:N, N+1] = sin_psi + M_t_vn @ Vn_lz
    rhs[0:N] = -xi * dfc_ds - M_t_vn @ Vn_rhs

    wrh = r * h * du
    M[N, 0:N] = wrh @ Vn_vs
    M[N, N] = float(np.dot(wrh, Vn_p))
    M[N, N+1] = float(np.dot(wrh, Vn_lz))
    rhs[N] = -float(np.dot(wrh, Vn_rhs))

    M[N+1, 0:N] = h * du
    M[N+1, N] = 0.0
    M[N+1, N+1] = 0.0
    rhs[N+1] = 0.0

    sol = np.linalg.solve(M, rhs)
    vs = sol[0:N]
    p = float(sol[N])
    lam_z = float(sol[N+1])
    vn = Vn_rhs + Vn_vs @ vs + Vn_p * p + Vn_lz * lam_z
    return vs, vn, p, lam_z


def solve_tube_pde(inputs, settings):
    N = int(round(settings.get("N", 33)))
    if N % 2 == 0:
        N += 1
    rtol = float(settings.get("rtol", 1e-4))
    atol = float(settings.get("atol", 1e-6))
    asym_pert = float(settings.get("asym_pert", 0.02))

    L0 = inputs["L_over_2pi_r0"] * 2.0 * math.pi
    xi = inputs["xi_over_gamma"]
    kappa = inputs["kappa_tilde"]
    D_t = inputs["D_tilde"]
    k_t = inputs["k_tilde"]
    eta_ratio = inputs["eta_s_over_eta_b"]
    m = inputs["m"]
    pert_amp = inputs["pert_amp"]
    t_end = inputs["t_end_tilde"]

    du, D1, D2 = spectral_periodic_ops(N, L0)
    u = np.arange(N) * du

    def rhs(t, y):
        h, psi, r, c = y.reshape(4, N)
        inv_h = 1.0 / h
        Css = inv_h * (D1 @ psi)
        Cphi = np.sin(psi) / r
        Gamma_u = h * np.cos(psi) / r
        vs, vn, p, lz = solve_hydro_tube(h, r, Css, Cphi, Gamma_u, c, xi, kappa, eta_ratio, m, D1, D2, du)
        dh_du = D1 @ h
        dvn_du = D1 @ vn
        dh = h * Css * vn
        dpsi = -inv_h * dvn_du
        dr = vn * np.sin(psi)
        Gamma_s = Gamma_u * inv_h
        d_cvs_ds = inv_h * (D1 @ (c * vs))
        dc_du = D1 @ c
        d2c_ds2 = (inv_h * inv_h) * (D2 @ c) - (dh_du * (inv_h ** 3)) * dc_du
        dc = -d_cvs_ds - Gamma_s * c * vs - c * vn * (Css + Cphi) + D_t * (d2c_ds2 + Gamma_s * inv_h * dc_du) - k_t * (c - 1.0)
        return np.concatenate([dh, dpsi, dr, dc])

    # Stop early if finite-time pinch-off approaches r -> 0 (e.g., kappa = 0 Fig. 3B)
    def pinch_event(t, y):
        r = y[2 * N : 3 * N]
        return float(np.min(r) - 0.05)

    pinch_event.terminal = True
    pinch_event.direction = -1

    c_init = 1.0 + pert_amp * np.cos(2.0 * math.pi * u / L0) + asym_pert * np.sin(4.0 * math.pi * u / L0)
    y0 = np.concatenate([np.ones(N), np.full(N, 0.5 * math.pi), np.ones(N), c_init])
    V0 = math.pi * L0

    sol = solve_ivp(rhs, (0.0, t_end), y0, method="BDF", rtol=rtol, atol=atol, events=pinch_event)
    if not sol.success and len(sol.t) < 2:
        raise RuntimeError(f"tube_active_surface_pde failed: {sol.message}")

    h, psi, r, c = sol.y[:, -1].reshape(4, N)
    Css = (D1 @ psi) / h
    Cphi = np.sin(psi) / r
    Gamma_u = h * np.cos(psi) / r
    vs, vn, p, lz = solve_hydro_tube(h, r, Css, Cphi, Gamma_u, c, xi, kappa, eta_ratio, m, D1, D2, du)

    cos_psi = np.cos(psi)
    sin_psi = np.sin(psi)
    denom = float(np.dot(cos_psi, cos_psi))
    v0 = -float(np.dot(vn, cos_psi)) / denom if denom > 1e-12 else 0.0
    vs_cm = vs - v0 * sin_psi
    # Orient co-moving flow positively when wave moves in negative z (matching Fig. 3D sign convention)
    if v0 > 0.0:
        vs_cm = -vs_cm
    V_end = float(math.pi * np.sum((r * r) * sin_psi * h) * du)

    return {
        "r": r, "c": c, "Css": Css, "Cphi": Cphi, "vs_cm": vs_cm, "v0": v0,
        "V_ratio": V_end / V0, "t_final": float(sol.t[-1]),
    }


def run(inputs, settings):
    res = solve_tube_pde(inputs, settings)
    return {
        "r_min_dimless": float(np.min(res["r"])),
        "r_max_dimless": float(np.max(res["r"])),
        "c_min_dimless": float(np.min(res["c"])),
        "c_max_dimless": float(np.max(res["c"])),
        "wave_speed_dimless": float(abs(res["v0"])),
        "vs_comov_min_dimless": float(np.min(res["vs_cm"])),
        "vs_comov_max_dimless": float(np.max(res["vs_cm"])),
        "Css_min_dimless": float(np.min(res["Css"])),
        "Css_max_dimless": float(np.max(res["Css"])),
        "Cphi_max_dimless": float(np.max(res["Cphi"])),
        "volume_ratio": float(res["V_ratio"]),
    }

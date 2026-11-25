import math
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import CubicSpline


def fd4_matrices_sphere(N, du):
    D1e = np.zeros((N, N))
    D1o = np.zeros((N, N))
    D2e = np.zeros((N, N))
    D2o = np.zeros((N, N))
    w1 = (1.0 / 12.0, -2.0 / 3.0, 0.0, 2.0 / 3.0, -1.0 / 12.0)
    w2 = (-1.0 / 12.0, 4.0 / 3.0, -2.5, 4.0 / 3.0, -1.0 / 12.0)
    offsets = (-2, -1, 0, 1, 2)
    for i in range(N):
        for off, c1, c2 in zip(offsets, w1, w2):
            j = i + off
            sgn_odd = 1.0
            if j < 0:
                j = -j
                sgn_odd = -1.0
            elif j >= N:
                j = 2 * (N - 1) - j
                sgn_odd = -1.0
            if c1 != 0.0:
                D1e[i, j] += c1 / du
                D1o[i, j] += (c1 / du) * sgn_odd
            if c2 != 0.0:
                D2e[i, j] += c2 / (du * du)
                D2o[i, j] += (c2 / (du * du)) * sgn_odd
    return D1e, D1o, D2e, D2o


def simpson_weights(N, du):
    w = np.ones(N)
    if N % 2 == 1 and N >= 3:
        w[1:-1:2] = 4.0
        w[2:-2:2] = 2.0
        w *= du / 3.0
    else:
        w[0] = 0.5
        w[-1] = 0.5
        w *= du
    return w


def init_spheroid(e, v_ratio, N):
    scale = (3.0 * v_ratio / (4.0 * math.pi)) ** (1.0 / 3.0)
    if e == 0.0:
        L0 = math.pi * scale
        u = np.linspace(0.0, L0, N)
        th = u / scale
        h = np.ones(N)
        r = scale * np.sin(th)
        r[0] = 0.0
        r[-1] = 0.0
        Css = np.full(N, 1.0 / scale)
        Cphi = np.full(N, 1.0 / scale)
        Gamma_u = np.zeros(N)
        Gamma_u[1:-1] = np.cos(th[1:-1]) / r[1:-1]
        return L0, u, h, r, Css, Cphi, Gamma_u

    b = scale * ((1.0 - e * e) ** (-1.0 / 3.0))
    a = scale * ((1.0 - e * e) ** (1.0 / 6.0))
    th_dense = np.linspace(0.0, math.pi, 4001)
    dr_d = a * np.cos(th_dense)
    dz_d = b * np.sin(th_dense)
    ds_dth = np.sqrt(dr_d * dr_d + dz_d * dz_d)
    # High-order cumulative integral via CubicSpline antiderivative
    s_spline = CubicSpline(th_dense, ds_dth).antiderivative()
    s_d = s_spline(th_dense) - s_spline(0.0)
    L0 = float(s_d[-1])
    u = np.linspace(0.0, L0, N)
    th = CubicSpline(s_d, th_dense)(u)
    th[0] = 0.0
    th[-1] = math.pi
    r = a * np.sin(th)
    r[0] = 0.0
    r[-1] = 0.0
    dr = a * np.cos(th)
    dz = b * np.sin(th)
    ds = np.sqrt(dr * dr + dz * dz)
    cos_psi = dr / ds
    Css = (a * b) / (ds ** 3)
    Cphi = b / (a * ds)
    Gamma_u = np.zeros(N)
    Gamma_u[1:-1] = cos_psi[1:-1] / r[1:-1]
    h = np.ones(N)
    return L0, u, h, r, Css, Cphi, Gamma_u


def solve_hydro_sphere(h, r, Css, Cphi, Gamma_u, c, xi, kappa, eta_ratio, m, D1e, D1o, D2e, D2o, w_simp, pole_kappa_coeff=8.0):
    """Solve axisymmetric spherical Eqs. (S32), (S33), (S41), (S54), (S59) with vn[0]=vn[-1].
    pole_kappa_coeff = 8.0 matches printed Eq. (S41); 16.0/3.0 is the exact pole Taylor limit.
    """
    N = len(h)
    ep = eta_ratio + 1.0
    em = 1.0 - eta_ratio
    cm = c ** m
    fc = cm / (1.0 + cm)
    Ckk = Css + Cphi
    inv_h = 1.0 / h
    Gamma_s = Gamma_u * inv_h
    sin_psi = r * Cphi
    cos_psi = r * Gamma_s
    cos_psi[0] = 1.0
    cos_psi[-1] = -1.0

    Ds_e = inv_h[:, None] * D1e
    Ds_o = inv_h[:, None] * D1o
    dh_du = D1e @ h
    dCss_ds = Ds_e @ Css
    dCphi_ds = np.zeros(N)
    dCphi_ds[1:-1] = Gamma_s[1:-1] * (Css[1:-1] - Cphi[1:-1])
    dGamma_ds = np.zeros(N)
    dGamma_ds[1:-1] = -Cphi[1:-1] * Css[1:-1] - Gamma_s[1:-1] ** 2
    dfc_ds = Ds_e @ fc

    # Eliminate vn algebraically from Eq. (S33) & (S41)
    A_n = ep * (Css * Css + Cphi * Cphi) + 2.0 * em * Css * Cphi
    A_n[0] = 4.0 * Css[0] * Css[0]
    A_n[-1] = 4.0 * Css[-1] * Css[-1]
    inv_An = 1.0 / A_n

    B_n = (ep * Css + em * Cphi)[:, None] * Ds_o
    diag_idx = np.arange(N)
    B_n[diag_idx, diag_idx] += (ep * Cphi + em * Css) * Gamma_s
    B_n[0, :] = 4.0 * Css[0] * Ds_o[0, :]
    B_n[-1, :] = 4.0 * Css[-1] * Ds_o[-1, :]

    rhs_n = -(xi * fc + 1.0) * Ckk
    if kappa != 0.0:
        Ds2_e = (inv_h * inv_h)[:, None] * D2e - (dh_du * (inv_h ** 3))[:, None] * D1e
        rhs_n = rhs_n + kappa * Ckk * ((Css - Cphi) ** 2) + 2.0 * kappa * (Ds2_e @ Ckk + Gamma_s * (Ds_e @ Ckk))
        d2Css = Ds2_e @ Css
        rhs_n[0] = -2.0 * (xi * fc[0] + 1.0) * Css[0] + pole_kappa_coeff * kappa * d2Css[0]
        rhs_n[-1] = -2.0 * (xi * fc[-1] + 1.0) * Css[-1] + pole_kappa_coeff * kappa * d2Css[-1]
    else:
        rhs_n[0] = -2.0 * (xi * fc[0] + 1.0) * Css[0]
        rhs_n[-1] = -2.0 * (xi * fc[-1] + 1.0) * Css[-1]

    Vn_vs = -inv_An[:, None] * B_n
    Vn_p = inv_An
    Vn_lz = inv_An * cos_psi
    Vn_rhs = inv_An * rhs_n

    Ds_of_Dsvs = Ds_e @ Ds_o
    M_t_vs = ep * Ds_of_Dsvs + (em + 2.0 * eta_ratio) * Gamma_s[:, None] * Ds_o
    M_t_vs[diag_idx, diag_idx] += em * dGamma_ds - 2.0 * eta_ratio * (Gamma_s ** 2)

    M_t_vn = (ep * Css + em * Cphi)[:, None] * Ds_e
    M_t_vn[diag_idx, diag_idx] += ep * dCss_ds + em * dCphi_ds + 2.0 * eta_ratio * Gamma_s * (Css - Cphi)

    M = np.zeros((N + 2, N + 2))
    rhs = np.zeros(N + 2)
    M[1:N-1, 0:N] = M_t_vs[1:N-1, :] + M_t_vn[1:N-1, :] @ Vn_vs
    M[1:N-1, N] = M_t_vn[1:N-1, :] @ Vn_p
    M[1:N-1, N+1] = sin_psi[1:N-1] + M_t_vn[1:N-1, :] @ Vn_lz
    rhs[1:N-1] = -xi * dfc_ds[1:N-1] - M_t_vn[1:N-1, :] @ Vn_rhs

    M[0, 0] = 1.0
    M[N-1, N-1] = 1.0

    wrh = w_simp * r * h
    M[N, 0:N] = wrh @ Vn_vs
    M[N, N] = float(np.dot(wrh, Vn_p))
    M[N, N+1] = float(np.dot(wrh, Vn_lz))
    rhs[N] = -float(np.dot(wrh, Vn_rhs))

    M[N+1, 0:N] = Vn_vs[0, :] - Vn_vs[-1, :]
    M[N+1, N] = Vn_p[0] - Vn_p[-1]
    M[N+1, N+1] = Vn_lz[0] - Vn_lz[-1]
    rhs[N+1] = -(Vn_rhs[0] - Vn_rhs[-1])

    sol = np.linalg.solve(M, rhs)
    vs = sol[0:N]
    p = float(sol[N])
    lam_z = float(sol[N+1])
    vn = Vn_rhs + Vn_vs @ vs + Vn_p * p + Vn_lz * lam_z
    return vs, vn, p, lam_z


def solve_single_step(inputs, settings, pole_kappa_coeff=8.0):
    N = int(round(settings.get("N", 101)))
    if N % 2 == 0:
        N += 1
    L0, u, h, r, Css, Cphi, Gamma_u = init_spheroid(inputs["e"], inputs["V_over_R03"], N)
    du = L0 / (N - 1)
    D1e, D1o, D2e, D2o = fd4_matrices_sphere(N, du)
    w_simp = simpson_weights(N, du)
    c = np.ones(N)
    vs, vn, p, _ = solve_hydro_sphere(
        h, r, Css, Cphi, Gamma_u, c,
        inputs["xi_over_gamma"], inputs["kappa_tilde"], inputs["eta_s_over_eta_b"], inputs["m"],
        D1e, D1o, D2e, D2o, w_simp, pole_kappa_coeff=pole_kappa_coeff,
    )
    return vs, vn, p


def run(inputs, settings):
    vs, vn, p = solve_single_step(inputs, settings, pole_kappa_coeff=8.0)
    return {
        "vn_max_dimless": float(np.max(np.abs(vn))),
        "vs_max_dimless": float(np.max(np.abs(vs))),
        "vn_pole_dimless": float(vn[0]),
        "p_dimless": float(p),
    }

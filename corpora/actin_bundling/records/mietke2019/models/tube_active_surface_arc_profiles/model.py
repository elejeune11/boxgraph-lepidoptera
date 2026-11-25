import importlib.util
import math
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

_base = Path(__file__).resolve().parent.parent
_spec = importlib.util.spec_from_file_location("_tube_pde_prof", _base / "tube_active_surface_pde" / "model.py")
_tube = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_tube)


def _periodic_spline(s_norm, y):
    s_ext = np.concatenate([s_norm, [1.0]])
    y_ext = np.concatenate([y, [y[0]]])
    return CubicSpline(s_ext, y_ext, bc_type="periodic")


def run(inputs, settings):
    N = int(round(settings.get("N", 33)))
    if N % 2 == 0:
        N += 1
    rtol = float(settings.get("rtol", 1e-4))
    atol = float(settings.get("atol", 1e-6))
    asym_pert = float(settings.get("asym_pert", 0.02))
    max_calls = int(settings.get("max_calls", 6000))

    L_rat = float(inputs["L_over_2pi_r0"])
    xi = float(inputs["xi_over_gamma"])
    kappa = float(inputs["kappa_tilde"])
    D_t = float(inputs["D_tilde"])
    k_t = float(inputs["k_tilde"])
    eta_ratio = float(inputs["eta_s_over_eta_b"])
    m = float(inputs["m"])
    pert_amp = float(inputs["pert_amp"])
    t_end = float(inputs["t_end_tilde"])
    c_peak_target = float(inputs["c_peak_s_over_L"])
    xs = [float(x) for x in inputs["sample_curve"]["x"]]

    if L_rat <= 0.0 or kappa < 0.0 or D_t < 0.0 or k_t < 0.0 or xi < 0.0 or eta_ratio <= 0.0 or m <= 0.0 or t_end < 0.0:
        raise ValueError("Unphysical negative/non-positive tubular PDE parameter")

    L0 = L_rat * 2.0 * math.pi
    du, D1, D2 = _tube.spectral_periodic_ops(N, L0)
    u = np.arange(N) * du
    D1T = D1.T
    D2T = D2.T
    ep = eta_ratio + 1.0
    em = 1.0 - eta_ratio
    diag_idx = np.arange(N)
    calls = [0]

    def rhs_vec(t, Y):
        calls[0] += 1
        if calls[0] > max_calls:
            raise RuntimeError(f"Tubular BDF evaluation budget ({max_calls}) exceeded")
        single = Y.ndim == 1
        if single:
            Y = Y[:, None]
        YT = Y.T.reshape(-1, 4, N)
        K = YT.shape[0]
        h = YT[:, 0, :]
        psi = YT[:, 1, :]
        r = YT[:, 2, :]
        c = YT[:, 3, :]
        if np.any(h <= 0.05) or np.any(r <= 0.03) or np.any(c < -0.5) or not np.all(np.isfinite(YT)):
            raise RuntimeError("Tubular solution left physical domain")

        inv_h = 1.0 / h
        Css = inv_h * (psi @ D1T)
        sin_psi = np.sin(psi)
        cos_psi = np.cos(psi)
        Cphi = sin_psi / r
        Gamma_s = cos_psi / r
        Ckk = Css + Cphi

        Ds = inv_h[:, :, None] * D1[None, :, :]
        dh_du = h @ D1T
        dCss_ds = inv_h * (Css @ D1T)
        dCphi_ds = Gamma_s * (Css - Cphi)
        dGamma_ds = -Cphi * Css - Gamma_s * Gamma_s

        if xi != 0.0:
            c_pos = np.maximum(c, 0.0)
            cm = c_pos ** m
            fc = cm / (1.0 + cm)
            dfc_ds = inv_h * (fc @ D1T)
        else:
            fc = np.zeros_like(c)
            dfc_ds = np.zeros_like(c)

        A_n = ep * (Css * Css + Cphi * Cphi) + 2.0 * em * Css * Cphi
        inv_An = 1.0 / A_n
        B_n = (ep * Css + em * Cphi)[:, :, None] * Ds
        B_n[:, diag_idx, diag_idx] += (ep * Cphi + em * Css) * Gamma_s

        rhs_n = -(xi * fc + 1.0) * Ckk
        if kappa != 0.0:
            dCkk_du = Ckk @ D1T
            d2Ckk = (inv_h * inv_h) * (Ckk @ D2T) - (dh_du * (inv_h ** 3)) * dCkk_du
            rhs_n = rhs_n + kappa * Ckk * ((Css - Cphi) ** 2) + 2.0 * kappa * (d2Ckk + Gamma_s * inv_h * dCkk_du)

        Vn_vs = -inv_An[:, :, None] * B_n
        Vn_p = inv_An
        Vn_lz = inv_An * cos_psi
        Vn_rhs = inv_An * rhs_n

        Ds2 = (inv_h * inv_h)[:, :, None] * D2[None, :, :] - (dh_du * (inv_h ** 3))[:, :, None] * D1[None, :, :]
        M_t_vs = ep * Ds2 + (em + 2.0 * eta_ratio) * Gamma_s[:, :, None] * Ds
        M_t_vs[:, diag_idx, diag_idx] += em * dGamma_ds - 2.0 * eta_ratio * (Gamma_s ** 2)

        M_t_vn = (ep * Css + em * Cphi)[:, :, None] * Ds
        M_t_vn[:, diag_idx, diag_idx] += ep * dCss_ds + em * dCphi_ds + 2.0 * eta_ratio * Gamma_s * (Css - Cphi)

        M = np.zeros((K, N + 2, N + 2))
        rhs_sys = np.zeros((K, N + 2))
        M[:, 0:N, 0:N] = M_t_vs + np.matmul(M_t_vn, Vn_vs)
        M[:, 0:N, N] = np.matmul(M_t_vn, Vn_p[:, :, None])[:, :, 0]
        M[:, 0:N, N + 1] = sin_psi + np.matmul(M_t_vn, Vn_lz[:, :, None])[:, :, 0]
        rhs_sys[:, 0:N] = -xi * dfc_ds - np.matmul(M_t_vn, Vn_rhs[:, :, None])[:, :, 0]

        wrh = r * h * du
        M[:, N, 0:N] = np.matmul(wrh[:, None, :], Vn_vs)[:, 0, :]
        M[:, N, N] = np.sum(wrh * Vn_p, axis=1)
        M[:, N, N + 1] = np.sum(wrh * Vn_lz, axis=1)
        rhs_sys[:, N] = -np.sum(wrh * Vn_rhs, axis=1)

        M[:, N + 1, 0:N] = h * du

        sol = np.linalg.solve(M, rhs_sys[:, :, None])[:, :, 0]
        vs = sol[:, 0:N]
        p = sol[:, N]
        lam_z = sol[:, N + 1]
        vn = Vn_rhs + np.matmul(Vn_vs, vs[:, :, None])[:, :, 0] + Vn_p * p[:, None] + Vn_lz * lam_z[:, None]

        dvn_du = vn @ D1T
        dh = h * Css * vn
        dpsi = -inv_h * dvn_du
        dr = vn * sin_psi
        d_cvs_ds = inv_h * ((c * vs) @ D1T)
        dc_du = c @ D1T
        d2c_ds2 = (inv_h * inv_h) * (c @ D2T) - (dh_du * (inv_h ** 3)) * dc_du
        dc = -d_cvs_ds - Gamma_s * c * vs - c * vn * Ckk + D_t * (d2c_ds2 + Gamma_s * inv_h * dc_du) - k_t * (c - 1.0)

        out = np.concatenate([dh, dpsi, dr, dc], axis=1).T
        return out[:, 0] if single else out

    c_init = 1.0 + pert_amp * np.cos(2.0 * math.pi * u / L0) + asym_pert * np.sin(4.0 * math.pi * u / L0)
    y0 = np.concatenate([np.ones(N), np.full(N, 0.5 * math.pi), np.ones(N), c_init])
    sol = solve_ivp(rhs_vec, (0.0, t_end), y0, method="BDF", rtol=rtol, atol=atol, vectorized=True)
    if not sol.success:
        raise RuntimeError(f"tube_active_surface_arc_profiles failed: {sol.message}")

    h, psi, r, c = sol.y[:, -1].reshape(4, N)
    Css = (D1 @ psi) / h
    Cphi = np.sin(psi) / r
    Gamma_u = h * np.cos(psi) / r
    vs, vn, p, lz = _tube.solve_hydro_tube(h, r, Css, Cphi, Gamma_u, c, xi, kappa, eta_ratio, m, D1, D2, du)

    cos_psi = np.cos(psi)
    sin_psi = np.sin(psi)
    denom = float(np.dot(cos_psi, cos_psi))
    v0 = -float(np.dot(vn, cos_psi)) / denom if denom > 1e-12 else 0.0
    vs_cm = vs - v0 * sin_psi

    u_ext = np.concatenate([u, [L0]])
    h_ext = np.concatenate([h, [h[0]]])
    s_ant = CubicSpline(u_ext, h_ext, bc_type="periodic").antiderivative()
    s_ext = s_ant(u_ext) - s_ant(0.0)
    s_norm = s_ext[:-1] / s_ext[-1]

    if v0 > 0.0:
        s_norm = np.concatenate([[0.0], (1.0 - s_norm[1:])[::-1]])
        c = np.concatenate([[c[0]], c[1:][::-1]])
        vs_cm = np.concatenate([[-vs_cm[0]], (-vs_cm[1:])[::-1]])
        Css = np.concatenate([[Css[0]], Css[1:][::-1]])
        Cphi = np.concatenate([[Cphi[0]], Cphi[1:][::-1]])

    sp_c = _periodic_spline(s_norm, c)
    sp_vs = _periodic_spline(s_norm, vs_cm)
    sp_Css = _periodic_spline(s_norm, Css)
    sp_Cphi = _periodic_spline(s_norm, Cphi)

    grid_fine = np.linspace(0.0, 1.0, 4001, endpoint=False)
    s_peak = float(grid_fine[int(np.argmax(sp_c(grid_fine)))])
    shift = c_peak_target - s_peak
    x_eval = (np.array(xs, dtype=float) - shift) % 1.0

    return {
        "c_profile": {"x": xs, "y": [float(v) for v in sp_c(x_eval)]},
        "vs_profile": {"x": xs, "y": [float(v) for v in sp_vs(x_eval)]},
        "Css_profile": {"x": xs, "y": [float(v) for v in sp_Css(x_eval)]},
        "Cphi_profile": {"x": xs, "y": [float(v) for v in sp_Cphi(x_eval)]},
    }

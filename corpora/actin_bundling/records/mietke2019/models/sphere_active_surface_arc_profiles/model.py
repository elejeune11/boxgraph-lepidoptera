import importlib.util
import math
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

_base = Path(__file__).resolve().parent.parent

_hyd_spec = importlib.util.spec_from_file_location("_sphere_hyd_prof", _base / "sphere_hydrodynamics_single_step_eqs32_33" / "model.py")
_hyd = importlib.util.module_from_spec(_hyd_spec)
_hyd_spec.loader.exec_module(_hyd)


def _cumint4(y, u):
    sp = CubicSpline(u, y).antiderivative()
    return sp(u) - sp(u[0])


def solve_with_c_init(inputs, settings):
    N = int(round(settings.get("N", 61)))
    if N % 2 == 0:
        N += 1
    rtol = float(settings.get("rtol", 1e-5))
    atol = float(settings.get("atol", 1e-7))
    max_nfev = int(settings.get("max_nfev", 6000))

    e_init = float(inputs["e_init"])
    v_ratio = float(inputs["V_over_R03"])
    xi = float(inputs["xi_over_gamma"])
    kappa = float(inputs["kappa_tilde"])
    D_t = float(inputs["D_tilde"])
    k_t = float(inputs["k_tilde"])
    eta_ratio = float(inputs["eta_s_over_eta_b"])
    m = float(inputs["m"])
    pert_amp = float(inputs["pert_amp"])
    t_end = float(inputs["t_end_tilde"])
    c_init_mean = float(inputs.get("c_init_dimless", 1.0))

    if e_init < 0.0 or e_init >= 1.0:
        raise ValueError(f"Unphysical eccentricity e_init={e_init}")
    if kappa < 0.0 or D_t < 0.0 or k_t < 0.0 or xi < 0.0 or eta_ratio <= 0.0 or m <= 0.0 or v_ratio <= 0.0:
        raise ValueError("Unphysical negative/non-positive physical parameter")

    L0, u, h0, r0, Css0, Cphi0, Gu0 = _hyd.init_spheroid(e_init, v_ratio, N)
    du = L0 / (N - 1)
    D1e, D1o, D2e, D2o = _hyd.fd4_matrices_sphere(N, du)
    w_simp = _hyd.simpson_weights(N, du)
    c0_init = c_init_mean + pert_amp * np.cos(math.pi * u / L0)

    if t_end <= 0.0:
        vs, vn, _, _ = _hyd.solve_hydro_sphere(
            h0, r0, Css0, Cphi0, Gu0, c0_init, xi, kappa, eta_ratio, m,
            D1e, D1o, D2e, D2o, w_simp, pole_kappa_coeff=8.0,
        )
        return {
            "u_norm": u / L0, "s_norm": u / L0,
            "c": c0_init, "vs": vs, "vn": vn, "Css": Css0, "Cphi": Cphi0, "h": h0,
        }

    calls = [0]

    def rhs(t, y):
        calls[0] += 1
        if calls[0] > max_nfev:
            raise RuntimeError(f"BDF evaluation budget ({max_nfev}) exceeded")
        h, Css, Cphi, Gu, r, c = y.reshape(6, N)
        if np.any(h <= 0.05) or np.any(c < -0.5) or not np.all(np.isfinite(y)):
            raise RuntimeError("Solution left physical domain")
        vs, vn, p, lz = _hyd.solve_hydro_sphere(
            h, r, Css, Cphi, Gu, c, xi, kappa, eta_ratio, m,
            D1e, D1o, D2e, D2o, w_simp, pole_kappa_coeff=8.0,
        )
        inv_h = 1.0 / h
        dh_du = D1e @ h
        dvn_du = D1e @ vn
        d2vn_ds2 = inv_h * (D1o @ (inv_h * dvn_du))
        dh = h * Css * vn
        dCss = -(Css * Css) * vn - d2vn_ds2
        dCphi = -(Cphi * Cphi) * vn - (inv_h * inv_h) * Gu * dvn_du
        dCphi[0] = dCss[0]
        dCphi[-1] = dCss[-1]
        dGu = D1e @ (Cphi * vn)
        dGu[0] = 0.0
        dGu[-1] = 0.0
        dr = vn * r * Cphi
        dr[0] = 0.0
        dr[-1] = 0.0
        Gamma_s = Gu * inv_h
        d_cvs_ds = inv_h * (D1o @ (c * vs))
        dc_du = D1e @ c
        dc_ds = inv_h * dc_du
        d2c_ds2 = (inv_h * inv_h) * (D2e @ c) - (dh_du * (inv_h ** 3)) * dc_du
        dc = -d_cvs_ds - Gamma_s * c * vs - c * vn * (Css + Cphi) + D_t * (d2c_ds2 + Gamma_s * dc_ds) - k_t * (c - 1.0)
        dc[0] = 2.0 * (-d_cvs_ds[0] - c[0] * vn[0] * Css[0] + D_t * d2c_ds2[0]) - k_t * (c[0] - 1.0)
        dc[-1] = 2.0 * (-d_cvs_ds[-1] - c[-1] * vn[-1] * Css[-1] + D_t * d2c_ds2[-1]) - k_t * (c[-1] - 1.0)
        return np.concatenate([dh, dCss, dCphi, dGu, dr, dc])

    y0 = np.concatenate([h0, Css0, Cphi0, Gu0, r0, c0_init])
    sol = solve_ivp(rhs, (0.0, t_end), y0, method="BDF", rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(f"sphere_active_surface_arc_profiles failed: {sol.message}")
    h, Css, Cphi, Gu, r, c = sol.y[:, -1].reshape(6, N)
    vs, vn, _, _ = _hyd.solve_hydro_sphere(
        h, r, Css, Cphi, Gu, c, xi, kappa, eta_ratio, m,
        D1e, D1o, D2e, D2o, w_simp, pole_kappa_coeff=8.0,
    )
    s_arc = _cumint4(h, u)
    return {
        "u_norm": u / L0, "s_norm": s_arc / s_arc[-1],
        "c": c, "vs": vs, "vn": vn, "Css": Css, "Cphi": Cphi, "h": h,
    }


def run(inputs, settings):
    xs = [float(x) for x in inputs["sample_curve"]["x"]]
    res = solve_with_c_init(inputs, settings)
    s_norm = res["s_norm"]
    x_arr = np.array(xs, dtype=float)
    return {
        "c_profile": {"x": xs, "y": [float(v) for v in CubicSpline(s_norm, res["c"])(x_arr)]},
        "vs_profile": {"x": xs, "y": [float(v) for v in CubicSpline(s_norm, res["vs"])(x_arr)]},
        "vn_profile": {"x": xs, "y": [float(v) for v in CubicSpline(s_norm, res["vn"])(x_arr)]},
        "Css_profile": {"x": xs, "y": [float(v) for v in CubicSpline(s_norm, res["Css"])(x_arr)]},
        "Cphi_profile": {"x": xs, "y": [float(v) for v in CubicSpline(s_norm, res["Cphi"])(x_arr)]},
    }

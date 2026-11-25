import importlib.util
import math
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

_hydro_path = Path(__file__).resolve().parent.parent / "sphere_hydrodynamics_single_step_eqs32_33" / "model.py"
_spec = importlib.util.spec_from_file_location("_sphere_hydro_pde", _hydro_path)
_hydro = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_hydro)


def _cumint4(y, u):
    sp = CubicSpline(u, y).antiderivative()
    return sp(u) - sp(u[0])


def solve_sphere_pde(inputs, settings, pole_kappa_coeff=8.0):
    N = int(round(settings.get("N", 61)))
    if N % 2 == 0:
        N += 1
    rtol = float(settings.get("rtol", 1e-5))
    atol = float(settings.get("atol", 1e-7))

    e_init = inputs["e_init"]
    v_ratio = inputs["V_over_R03"]
    xi = inputs["xi_over_gamma"]
    kappa = inputs["kappa_tilde"]
    D_t = inputs["D_tilde"]
    k_t = inputs["k_tilde"]
    eta_ratio = inputs["eta_s_over_eta_b"]
    m = inputs["m"]
    pert_amp = inputs["pert_amp"]
    t_end = inputs["t_end_tilde"]

    L0, u, h0, r0, Css0, Cphi0, Gu0 = _hydro.init_spheroid(e_init, v_ratio, N)
    du = L0 / (N - 1)
    D1e, D1o, D2e, D2o = _hydro.fd4_matrices_sphere(N, du)
    w_simp = _hydro.simpson_weights(N, du)

    # Polar l=1 perturbation cos(pi*u/L0) has zero area-weighted mean on sphere
    c0_init = 1.0 + pert_amp * np.cos(math.pi * u / L0)
    psi0 = _cumint4(Css0 * h0, u)
    V0 = float(math.pi * np.sum(w_simp * (r0 ** 2) * np.sin(psi0) * h0))
    Nc0 = float(2.0 * math.pi * np.sum(w_simp * c0_init * r0 * h0))

    def rhs(t, y):
        h, Css, Cphi, Gu, r, c = y.reshape(6, N)
        vs, vn, p, lz = _hydro.solve_hydro_sphere(
            h, r, Css, Cphi, Gu, c, xi, kappa, eta_ratio, m,
            D1e, D1o, D2e, D2o, w_simp, pole_kappa_coeff=pole_kappa_coeff,
        )
        inv_h = 1.0 / h
        dh_du = D1e @ h
        dvn_du = D1e @ vn
        # Exact Eq. (24) / Eq. (S49): -(1/h) d_u((1/h) d_u vn)
        d2vn_ds2 = inv_h * (D1o @ (inv_h * dvn_du))

        # Eqs. (S46), (S49), (S50), (S51), (S53)
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

        # Eqs. (S34) and (S42)
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
        raise RuntimeError(f"sphere_active_surface_pde failed: {sol.message}")

    h, Css, Cphi, Gu, r, c = sol.y[:, -1].reshape(6, N)
    vs, vn, p, lz = _hydro.solve_hydro_sphere(
        h, r, Css, Cphi, Gu, c, xi, kappa, eta_ratio, m,
        D1e, D1o, D2e, D2o, w_simp, pole_kappa_coeff=pole_kappa_coeff,
    )
    # Reconstruct arc length and tangent angle Eqs. (17) and check Eqs. (S57), (S62)
    s_arc = _cumint4(h, u)
    psi_rec = _cumint4(Css * h, u)
    V_end = float(math.pi * np.sum(w_simp * (r ** 2) * np.sin(psi_rec) * h))
    Nc_end = float(2.0 * math.pi * np.sum(w_simp * c * r * h))
    s_norm = s_arc / s_arc[-1]

    return {
        "s_norm": s_norm, "c": c, "vs": vs, "vn": vn, "Css": Css, "Cphi": Cphi, "h": h, "r": r,
        "V_ratio": V_end / V0, "Nc_ratio": Nc_end / Nc0,
    }


def run(inputs, settings):
    res = solve_sphere_pde(inputs, settings, pole_kappa_coeff=8.0)
    c = res["c"]
    vs = res["vs"]
    vn = res["vn"]
    Css = res["Css"]
    Cphi = res["Cphi"]
    return {
        "c_pole_dimless": float(c[0]),
        "c_min_dimless": float(np.min(c)),
        "c_max_dimless": float(np.max(c)),
        "vs_min_dimless": float(np.min(vs)),
        "Css_pole_dimless": float(Css[0]),
        "Css_max_dimless": float(np.max(Css)),
        "Cphi_max_dimless": float(np.max(Cphi)),
        "Cphi_opposite_pole_dimless": float(Cphi[-1]),
        "vn_max_dimless": float(np.max(np.abs(vn))),
        "volume_ratio": float(res["V_ratio"]),
        "molecule_ratio": float(res["Nc_ratio"]),
        "rel_vol_error": float(abs(res["V_ratio"] - 1.0)),
        "rel_mass_error": float(abs(res["Nc_ratio"] - 1.0)),
    }

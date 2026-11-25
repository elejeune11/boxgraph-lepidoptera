import importlib.util
import math
from pathlib import Path
import numpy as np

_tube_path = Path(__file__).resolve().parent.parent / "tube_active_surface_pde" / "model.py"
_spec = importlib.util.spec_from_file_location("_tube_pde_s2b", _tube_path)
_tube = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_tube)


def _evolve_euler(L0, xi, kappa, eta_ratio, t_end, N, Nt):
    du, D1, D2 = _tube.spectral_periodic_ops(N, L0)
    u = np.arange(N) * du
    k = 2.0 * math.pi / L0
    eps = 0.05
    # Smooth initial perturbed tube r(u) = 1 + eps*cos(k*u), z(u) = u
    r = 1.0 + eps * np.cos(k * u)
    dr_du = -eps * k * np.sin(k * u)
    h = np.sqrt(1.0 + dr_du * dr_du)
    psi = np.arctan2(1.0, dr_du)
    c = np.ones(N)
    dt = t_end / float(Nt)

    for _ in range(Nt):
        inv_h = 1.0 / h
        Css = inv_h * (D1 @ psi)
        Cphi = np.sin(psi) / r
        Gamma_u = h * np.cos(psi) / r
        vs, vn, p, lz = _tube.solve_hydro_tube(h, r, Css, Cphi, Gamma_u, c, xi, kappa, eta_ratio, 2, D1, D2, du)
        dvn_du = D1 @ vn
        h = h + dt * (h * Css * vn)
        psi = psi + dt * (-inv_h * dvn_du)
        r = r + dt * (vn * np.sin(psi))

    Css = (D1 @ psi) / h
    Cphi = np.sin(psi) / r
    return Css, Cphi


def run(inputs, settings):
    N = int(round(settings.get("N", 25)))
    if N % 2 == 0:
        N += 1
    Nt_base = int(round(settings.get("Nt_base", 50)))
    Nt_ref = int(round(settings.get("Nt_ref", 1600)))

    L0 = inputs["L_over_2pi_r0"] * 2.0 * math.pi
    xi = inputs["xi_over_gamma"]
    kappa = inputs["kappa_tilde"]
    eta_ratio = inputs["eta_s_over_eta_b"]
    t_end = inputs["t_end_tilde"]

    Css_ref, Cphi_ref = _evolve_euler(L0, xi, kappa, eta_ratio, t_end, N, Nt_ref)
    nts = [Nt_base, 2 * Nt_base, 4 * Nt_base]
    err_css, err_cphi = [], []
    for nt in nts:
        Css, Cphi = _evolve_euler(L0, xi, kappa, eta_ratio, t_end, N, nt)
        err_css.append(float(np.max(np.abs(Css - Css_ref))))
        err_cphi.append(float(np.max(np.abs(Cphi - Cphi_ref))))

    logNt = np.log(np.array(nts, dtype=float))
    p_css = float(-np.polyfit(logNt, np.log(np.maximum(err_css, 1e-16)), 1)[0])
    p_cphi = float(-np.polyfit(logNt, np.log(np.maximum(err_cphi, 1e-16)), 1)[0])
    return {"order_Css": p_css, "order_Cphi": p_cphi}

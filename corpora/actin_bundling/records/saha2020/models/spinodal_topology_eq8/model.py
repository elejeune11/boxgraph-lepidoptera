import math
import numpy as np
from scipy.ndimage import label
from scipy.optimize import root


def _d_entries(p1, p2, c11, c12, c21, c22, chi, chip, alpha):
    d11 = -(2.0 * (p1 - c11) ** 2 + 8.0 * (p1 - c11) * (p1 - c12) + 2.0 * (p1 - c12) ** 2 + 2.0 * p2 * p2 * chip)
    d22 = -(2.0 * (p2 - c21) ** 2 + 8.0 * (p2 - c21) * (p2 - c22) + 2.0 * (p2 - c22) ** 2 + 2.0 * p1 * p1 * chip)
    d12 = -((chi + alpha) + 4.0 * p1 * p2 * chip)
    d21 = -((chi - alpha) + 4.0 * p1 * p2 * chip)
    return d11, d12, d21, d22


def run(inputs, settings):
    c11 = float(inputs["c11"])
    c12 = float(inputs["c12"])
    c21 = float(inputs["c21"])
    c22 = float(inputs["c22"])
    chi = float(inputs["chi"])
    chip = float(inputs["chip"])
    alpha = float(inputs["alpha"])

    n_grid = int(settings.get("n_grid", 201))
    phi_min = float(settings.get("phi_min", -0.8))
    phi_max = float(settings.get("phi_max", 1.2))

    grid = np.linspace(phi_min, phi_max, n_grid)
    P1, P2 = np.meshgrid(grid, grid)
    d11, d12, d21, d22 = _d_entries(P1, P2, c11, c12, c21, c22, chi, chip, alpha)
    tr = d11 + d22
    det = d11 * d22 - d12 * d21
    disc = tr * tr - 4.0 * det
    in_ci = disc < 0.0
    re_max = np.where(disc >= 0.0, 0.5 * (tr + np.sqrt(np.maximum(disc, 0.0))), 0.5 * tr)
    in_cr = re_max > 0.0

    _, n_comp = label(in_cr)
    _, n_osc = label(in_cr & in_ci)

    if abs(chip) < 1e-12:
        c1 = abs(c12 - c11) / 2.0
        c2 = abs(c22 - c21) / 2.0
        alpha_split = math.sqrt(chi * chi + (4.0 * max(c1 * c1, c2 * c2)) ** 2)
    else:
        # Saddle-point pinch-off of the bulk spinodal det(D/q^2) = 0, grad det = 0
        # Note det = d11*d22 - (chi + 4*chip*p1*p2)^2 + alpha^2 = alpha^2 - f(p1, p2),
        # so the saddle of det=0 is a saddle of f(p1, p2) = (chi + 4*chip*p1*p2)^2 - d11*d22!
        def f_val(p1, p2):
            a11, _, _, a22 = _d_entries(p1, p2, c11, c12, c21, c22, chi, chip, 0.0)
            return (chi + 4.0 * chip * p1 * p2) ** 2 - a11 * a22

        def grad_f(xy):
            p1, p2 = xy
            eps = 1e-7
            gx = (f_val(p1 + eps, p2) - f_val(p1 - eps, p2)) / (2.0 * eps)
            gy = (f_val(p1, p2 + eps) - f_val(p1, p2 - eps)) / (2.0 * eps)
            return [gx, gy]

        best_alpha = float("nan")
        for guess in ([0.22, 0.15], [0.26, 0.18], [0.15, 0.15], [0.20, 0.10]):
            sol = root(grad_f, guess, tol=1e-12)
            if not sol.success:
                continue
            p1s, p2s = sol.x
            fv = f_val(p1s, p2s)
            if fv <= 0.0:
                continue
            # Check Hessian signature (must be a saddle point, det(H) < 0)
            eps = 1e-5
            f0 = fv
            fxx = (f_val(p1s + eps, p2s) - 2.0 * f0 + f_val(p1s - eps, p2s)) / (eps * eps)
            fyy = (f_val(p1s, p2s + eps) - 2.0 * f0 + f_val(p1s, p2s - eps)) / (eps * eps)
            fxy = (
                f_val(p1s + eps, p2s + eps)
                - f_val(p1s + eps, p2s - eps)
                - f_val(p1s - eps, p2s + eps)
                + f_val(p1s - eps, p2s - eps)
            ) / (4.0 * eps * eps)
            if fxx * fyy - fxy * fxy < 0.0:
                best_alpha = math.sqrt(fv)
                break
        alpha_split = best_alpha

    return {
        "n_components": float(n_comp),
        "n_osc_components": float(n_osc),
        "alpha_split": float(alpha_split),
    }

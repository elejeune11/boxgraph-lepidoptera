import importlib.util
import math
from pathlib import Path
import numpy as np
from scipy.optimize import brentq

_jac_path = Path(__file__).resolve().parent.parent / "cylinder_jacobian_eqs80_82" / "model.py"
_spec = importlib.util.spec_from_file_location("_cyl_jac", _jac_path)
_jac = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_jac)


def run(inputs, settings):
    alpha = inputs["alpha"]
    xi_f0 = inputs["xi_f0_over_gamma"]
    kappa_tilde = inputs["kappa_tilde"]
    D_tilde = inputs["D_tilde"]
    k_turn = inputs["k_tilde"]
    eta_ratio = inputs["eta_s_over_eta_b"]

    n_scan = int(round(settings.get("n_scan", 20000)))
    k_max = float(settings.get("k_max", 20.0))

    def max_re_fn(kt):
        J11, J12, J21, J22 = _jac.jacobian_entries(kt, alpha, xi_f0, kappa_tilde, D_tilde, k_turn, eta_ratio)
        tr = J11 + J22
        det = J11 * J22 - J12 * J21
        disc = tr * tr - 4.0 * det
        return 0.5 * (tr + math.sqrt(disc)) if disc >= 0.0 else 0.5 * tr

    def disc_fn(kt):
        J11, J12, J21, J22 = _jac.jacobian_entries(kt, alpha, xi_f0, kappa_tilde, D_tilde, k_turn, eta_ratio)
        tr = J11 + J22
        det = J11 * J22 - J12 * J21
        return tr * tr - 4.0 * det

    ks = np.linspace(1e-4, k_max, n_scan)
    k2 = ks * ks
    bc = 4.0 * k2 + (2.0 * k2 - 3.0) * (k2 - 1.0)
    term = -kappa_tilde * bc + (1.0 + xi_f0) * (1.0 - k2)
    nu_bar = (eta_ratio + 1.0) / (4.0 * eta_ratio)
    J11 = nu_bar * term
    J12 = -0.5 * alpha
    J21 = -0.5 * term
    J22 = alpha - D_tilde * k2 - k_turn
    tr = J11 + J22
    det = J11 * J22 - J12 * J21
    disc = tr * tr - 4.0 * det
    unstable = (tr > 0.0) | (det < 0.0)

    if not np.any(unstable):
        return {
            "L_crit_over_2pi_r0": 0.0,
            "L_crit_over_r0": 0.0,
            "L_osc_min_over_2pi_r0": 0.0,
            "L_osc_max_over_2pi_r0": 0.0,
        }

    # Largest unstable wavenumber k_crit -> smallest unstable aspect ratio L_crit/(2*pi*r0) = 1/k_crit
    idx_last = int(np.where(unstable)[0][-1])
    if idx_last + 1 < len(ks):
        k_crit = brentq(max_re_fn, ks[idx_last], ks[idx_last + 1], xtol=1e-12)
    else:
        k_crit = float(ks[idx_last])
    L_crit = 1.0 / k_crit

    # Oscillatory unstable region: unstable and disc < 0
    osc = unstable & (disc < 0.0)
    if np.any(osc):
        osc_idx = np.where(osc)[0]
        i_lo = int(osc_idx[0])
        i_hi = int(osc_idx[-1])
        # Refine k_osc_lo (gives L_osc_max) and k_osc_hi (gives L_osc_min)
        if i_lo > 0 and disc[i_lo - 1] >= 0.0:
            k_osc_lo = brentq(disc_fn, ks[i_lo - 1], ks[i_lo], xtol=1e-12)
        else:
            k_osc_lo = float(ks[i_lo])
        if i_hi + 1 < len(ks) and not unstable[i_hi + 1]:
            k_osc_hi = k_crit
        elif i_hi + 1 < len(ks) and disc[i_hi + 1] >= 0.0:
            k_osc_hi = brentq(disc_fn, ks[i_hi], ks[i_hi + 1], xtol=1e-12)
        else:
            k_osc_hi = float(ks[i_hi])
        L_osc_min = 1.0 / k_osc_hi
        L_osc_max = 1.0 / k_osc_lo
    else:
        L_osc_min = 0.0
        L_osc_max = 0.0

    return {
        "L_crit_over_2pi_r0": L_crit,
        "L_crit_over_r0": 2.0 * math.pi * L_crit,
        "L_osc_min_over_2pi_r0": L_osc_min,
        "L_osc_max_over_2pi_r0": L_osc_max,
    }

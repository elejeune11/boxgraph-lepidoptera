import math
import numpy as np


def run(inputs, settings):
    aA = float(inputs["alpha_A"])
    aB = float(inputs["alpha_B"])
    k = float(inputs["kappa"])
    d = float(inputs["delta"])

    # Exact onset delta_c > kappa satisfies:
    # delta^2 * (delta^2 - kappa^2) = (0.5*(aA - aB)*delta + 0.5*k*(aA + aB))^2
    c_diff = 0.5 * (aA - aB)
    c_sum = 0.5 * k * (aA + aB)
    # Quartic: d^4 - (k^2 + c_diff^2)*d^2 - 2*c_diff*c_sum*d - c_sum^2 = 0
    poly = [1.0, 0.0, -(k * k + c_diff * c_diff), -2.0 * c_diff * c_sum, -c_sum * c_sum]
    roots = np.roots(poly)
    real_pos = [r.real for r in roots if abs(r.imag) < 1e-10 and r.real >= k - 1e-12]
    dc = min(real_pos) if real_pos else k
    dc_over_k = dc / k

    if d <= k or (aA + aB) >= 0.0:
        return {"delta_c_over_kappa": float(dc_over_k), "v": 0.0, "rho_A_t": None, "rho_B_t": None}

    shift = c_diff + c_sum / d
    v2 = d * d - k * k - shift * shift
    if v2 <= 0.0:
        return {"delta_c_over_kappa": float(dc_over_k), "v": 0.0, "rho_A_t": None, "rho_B_t": None}

    rho_A = math.sqrt(-(d - k) * (aA + aB) / (2.0 * d))
    rho_B = math.sqrt(-(d + k) * (aA + aB) / (2.0 * d))
    return {
        "delta_c_over_kappa": float(dc_over_k),
        "v": math.sqrt(v2),
        "rho_A_t": rho_A,
        "rho_B_t": rho_B,
    }


def symbolic():
    import sympy as sp

    alpha_A, alpha_B, kappa, delta = sp.symbols("alpha_A alpha_B kappa delta", real=True)
    shift = (alpha_A - alpha_B) / 2 + kappa * (alpha_A + alpha_B) / (2 * delta)
    # Symmetric case formula for delta_c_over_kappa when alpha_A = alpha_B = alpha:
    return {
        "v": sp.sqrt(delta**2 - kappa**2 - shift**2),
        "rho_A_t": sp.sqrt(-(delta - kappa) * (alpha_A + alpha_B) / (2 * delta)),
        "rho_B_t": sp.sqrt(-(delta + kappa) * (alpha_A + alpha_B) / (2 * delta)),
    }

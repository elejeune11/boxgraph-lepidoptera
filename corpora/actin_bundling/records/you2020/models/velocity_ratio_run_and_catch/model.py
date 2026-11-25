import math


def run(inputs, settings):
    aB = float(inputs["alpha_B"])
    k = float(inputs["kappa"])
    d = float(inputs["delta"])
    rA = float(inputs["rho_A"])
    rB = float(inputs["rho_B"])
    v_ratio_gen = -((k - d) * rB * rB) / ((k + d) * rA * rA)
    v_ratio_fs = (d * d - k * k) / (aB * aB)
    dc = math.sqrt(k * k + aB * aB)
    psi_rate = 2.0 * dc * (d - dc) / aB
    return {
        "v_ratio_general": v_ratio_gen,
        "v_ratio_fs": v_ratio_fs,
        "psi_growth_rate_near_dc": psi_rate,
        "sum_phase_slope_near_dc": 2.0 * aB,
    }


def symbolic():
    import sympy as sp

    alpha_B, kappa, delta, rho_A, rho_B = sp.symbols(
        "alpha_B kappa delta rho_A rho_B", real=True
    )
    dc = sp.sqrt(kappa**2 + alpha_B**2)
    return {
        "v_ratio_general": -((kappa - delta) * rho_B**2) / ((kappa + delta) * rho_A**2),
        "v_ratio_fs": (delta**2 - kappa**2) / alpha_B**2,
        "psi_growth_rate_near_dc": 2 * dc * (delta - dc) / alpha_B,
        "sum_phase_slope_near_dc": 2 * alpha_B,
    }

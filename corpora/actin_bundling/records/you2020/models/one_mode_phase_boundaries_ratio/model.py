import math


def run(inputs, settings):
    gA = float(inputs["gamma_A"])
    p0 = float(inputs["phi0_A"])
    aB = float(inputs["alpha_B"])
    k = float(inputs["kappa"])
    d = float(inputs["delta_over_kappa"]) * k
    abs_chi_diff = gA + p0 * p0 + (d * d - k * k) / aB
    abs_chi_osc = gA + p0 * p0 + aB
    dc_over_k = math.sqrt(1.0 + (aB / k) ** 2)
    return {
        "abs_chi_A_diff": abs_chi_diff,
        "abs_chi_A_osc": abs_chi_osc,
        "delta_c_over_kappa": dc_over_k,
    }


def symbolic():
    import sympy as sp

    gamma_A, phi0_A, alpha_B, kappa, delta_over_kappa = sp.symbols(
        "gamma_A phi0_A alpha_B kappa delta_over_kappa", real=True
    )
    delta = delta_over_kappa * kappa
    return {
        "abs_chi_A_diff": gamma_A + phi0_A**2 + (delta**2 - kappa**2) / alpha_B,
        "abs_chi_A_osc": gamma_A + phi0_A**2 + alpha_B,
        "delta_c_over_kappa": sp.sqrt(1 + alpha_B**2 / kappa**2),
    }

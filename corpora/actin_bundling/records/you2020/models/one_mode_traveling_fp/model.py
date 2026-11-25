import math


def run(inputs, settings):
    aA = float(inputs["alpha_A"])
    aB = float(inputs["alpha_B"])
    k = float(inputs["kappa"])
    d = float(inputs["delta"])
    dc = math.sqrt(k * k + aB * aB)
    dc_over_k = dc / k
    if d < dc - 1e-14 or (-aA - aB) < 0.0:
        return {
            "delta_c": dc,
            "delta_c_over_kappa": dc_over_k,
            "v": 0.0,
            "rho_A_t": None,
            "rho_B_t": None,
            "theta_t": None,
            "psi_t": None,
            "amp_ratio_t": None,
        }
    v = math.sqrt(max(0.0, d * d - dc * dc))
    rho_A = math.sqrt(-aA - aB)
    ratio = math.sqrt((d - k) / (d + k))
    rho_B = rho_A / ratio
    cos_arg = -min(1.0, math.sqrt(aB * aB / (d * d - k * k)))
    theta_t = math.acos(cos_arg)
    psi_t = math.pi - theta_t
    return {
        "delta_c": dc,
        "delta_c_over_kappa": dc_over_k,
        "v": v,
        "rho_A_t": rho_A,
        "rho_B_t": rho_B,
        "theta_t": theta_t,
        "psi_t": psi_t,
        "amp_ratio_t": ratio,
    }


def symbolic():
    import sympy as sp

    alpha_A, alpha_B, kappa, delta = sp.symbols("alpha_A alpha_B kappa delta", real=True)
    dc = sp.sqrt(kappa**2 + alpha_B**2)
    rho_A = sp.sqrt(-alpha_A - alpha_B)
    ratio = sp.sqrt((delta - kappa) / (delta + kappa))
    theta_t = sp.acos(-sp.sqrt(alpha_B**2 / (delta**2 - kappa**2)))
    return {
        "delta_c": dc,
        "delta_c_over_kappa": dc / kappa,
        "v": sp.sqrt(delta**2 - kappa**2 - alpha_B**2),
        "rho_A_t": rho_A,
        "rho_B_t": sp.sqrt((delta + kappa) / (delta - kappa)) * rho_A,
        "theta_t": theta_t,
        "psi_t": sp.pi - theta_t,
        "amp_ratio_t": ratio,
    }

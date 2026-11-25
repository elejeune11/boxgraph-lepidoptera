import math


def run(inputs, settings):
    abs_chi = float(inputs["abs_chi_A"])
    gA = float(inputs["gamma_A"])
    p0 = float(inputs["phi0_A"])
    aB = float(inputs["alpha_B"])
    k = float(inputs["kappa"])
    if abs_chi < gA + p0 * p0 + aB - 1e-12:
        return {"delta_c_over_kappa": None}
    return {"delta_c_over_kappa": math.sqrt(1.0 + (aB / k) ** 2)}

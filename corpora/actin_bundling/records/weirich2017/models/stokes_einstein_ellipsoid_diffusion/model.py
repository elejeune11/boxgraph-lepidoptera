import math

KB_J_PER_K = 1.380649e-23


def run(inputs, settings):
    L_m = inputs["L"] * 1e-9
    D_m2_s = KB_J_PER_K * inputs["T"] * math.log(2.0 * inputs["L"] / inputs["r"]) / (
        6.0 * math.pi * inputs["eta"] * L_m
    )
    return {"D": D_m2_s * 1e12}


def symbolic():
    import sympy as sp

    T, L, r, eta = sp.symbols("T L r eta", positive=True)
    # 1e12 (um^2/m^2) / 1e-9 (m/nm) = 1e21
    coeff = sp.Float(KB_J_PER_K) * sp.Integer(10**21)
    return {"D": coeff * T * sp.log(2 * L / r) / (6 * sp.pi * eta * L)}

import math

KB_J_PER_K = 1.380649e-23


def run(inputs, settings):
    L_m = inputs["L"] * 1e-9
    D_m2_s = inputs["D"] * 1e-12
    eta = KB_J_PER_K * inputs["T"] * math.log(2.0 * inputs["L"] / inputs["r"]) / (
        6.0 * math.pi * D_m2_s * L_m
    )
    return {"eta": eta}


def symbolic():
    import sympy as sp

    T, L, r, D = sp.symbols("T L r D", positive=True)
    coeff = sp.Float(KB_J_PER_K) * sp.Integer(10**21)
    return {"eta": coeff * T * sp.log(2 * L / r) / (6 * sp.pi * D * L)}

import math


def run(inputs, settings):
    V_bar = (inputs["V_s"] / inputs["psi"]) * (inputs["phibar"] - inputs["k_b"] / inputs["k"])
    if V_bar <= 0.0:
        return {"V_bar": None, "R_bar": None}
    R_bar = (3.0 * V_bar / (4.0 * math.pi)) ** (1.0 / 3.0)
    return {"V_bar": V_bar, "R_bar": R_bar}


def symbolic():
    import sympy as sp

    V_s, psi, phibar, k_b, k = sp.symbols("V_s psi phibar k_b k", positive=True)
    V_bar = (V_s / psi) * (phibar - k_b / k)
    R_bar = (3 * V_bar / (4 * sp.pi)) ** sp.Rational(1, 3)
    return {"V_bar": V_bar, "R_bar": R_bar}

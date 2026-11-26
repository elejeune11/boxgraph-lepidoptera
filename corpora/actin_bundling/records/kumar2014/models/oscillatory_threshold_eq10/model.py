import math


def run(inputs, settings):
    k1_l_sq = (inputs["omega"] * math.pi / inputs["L_over_l"]) ** 2
    num = (1.0 + inputs["alpha"]) * (1.0 + k1_l_sq * inputs["nu"])
    den = inputs["A0_fA"] + inputs["I0_fI"]
    return {"P_c": num / den}


def symbolic():
    import sympy as sp

    alpha, omega, L_over_l, nu = sp.symbols("alpha omega L_over_l nu", positive=True)
    A0_fA, I0_fI = sp.symbols("A0_fA I0_fI", real=True)
    return {"P_c": (1 + alpha) * (1 + omega**2 * sp.pi**2 * nu / L_over_l**2) / (A0_fA + I0_fI)}

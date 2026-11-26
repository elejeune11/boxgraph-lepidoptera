import math


def run(inputs, settings):
    return {
        "L_perp_max": inputs["a_perp"]
        * math.sqrt(inputs["N"])
        * math.exp(inputs["L"] * inputs["J"] * inputs["q"] / (4.0 * inputs["T"]))
    }


def symbolic():
    import sympy as sp
    a_perp, N, L, J, q, T = sp.symbols("a_perp N L J q T", positive=True)
    return {"L_perp_max": a_perp * sp.sqrt(N) * sp.exp(L * J * q / (4 * T))}

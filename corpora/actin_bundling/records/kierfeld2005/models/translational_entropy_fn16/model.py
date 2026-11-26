import math


def run(inputs, settings):
    N = inputs["N"]
    a_perp = inputs["a_perp"]
    L_perp = inputs["L_perp"]
    return {"S_tr": (N - 1.0) * math.log(N * (a_perp ** 2) / (L_perp ** 2))}


def symbolic():
    import sympy as sp
    N, a_perp, L_perp = sp.symbols("N a_perp L_perp", positive=True)
    return {"S_tr": (N - 1) * sp.log(N * a_perp ** 2 / L_perp ** 2)}

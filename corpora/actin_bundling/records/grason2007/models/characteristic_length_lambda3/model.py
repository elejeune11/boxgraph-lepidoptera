import math


def run(inputs, settings):
    return {"lambda_3perp": math.sqrt(inputs["K3_bar"] / inputs["mu_perp_bar"])}


def symbolic():
    import sympy as sp
    K3_bar, mu_perp_bar = sp.symbols("K3_bar mu_perp_bar", positive=True)
    return {"lambda_3perp": sp.sqrt(K3_bar / mu_perp_bar)}

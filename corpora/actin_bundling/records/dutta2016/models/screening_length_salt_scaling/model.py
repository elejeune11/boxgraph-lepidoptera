import math


def run(inputs, settings):
    # Unknown prefactor set to 1 (prefactor: unknown)
    return {"lambda_R_phys": 1.0 / math.sqrt(inputs["n"])}


def symbolic():
    import sympy as sp
    n = sp.Symbol("n", positive=True)
    return {"lambda_R_phys": 1 / sp.sqrt(n)}

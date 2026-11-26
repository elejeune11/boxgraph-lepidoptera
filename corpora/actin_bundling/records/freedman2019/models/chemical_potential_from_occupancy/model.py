import math


def run(inputs, settings):
    N = inputs["N"]
    mu = inputs["mu_kBT"]
    n_b = N / (1.0 + math.exp(-mu))
    mu_inv = -math.log(N / n_b - 1.0)
    return {"N_b": n_b, "mu_kBT_inv": mu_inv}


def symbolic():
    import sympy as sp

    N = sp.Symbol("N", positive=True)
    mu_kBT = sp.Symbol("mu_kBT", real=True)
    n_b = N / (1 + sp.exp(-mu_kBT))
    return {"N_b": n_b, "mu_kBT_inv": -sp.log(N / n_b - 1)}

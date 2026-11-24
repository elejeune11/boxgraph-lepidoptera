import math


def run(inputs, settings):
    k_b_fil = inputs["k_b_fil"]
    N = inputs["N"]
    chi = inputs["chi"]
    c_qj = inputs["c_qj"]
    alpha = inputs["alpha"]
    nl_term = (chi ** 2 * (N - 1.0)) / (1.0 + c_qj * (N + math.sqrt(N)) / alpha)
    return {"k_b": k_b_fil * N * (1.0 + nl_term)}


def symbolic():
    import sympy as sp

    k_b_fil, N, chi, c_qj, alpha = sp.symbols("k_b_fil N chi c_qj alpha", positive=True)
    return {"k_b": k_b_fil * N * (1 + chi ** 2 * (N - 1) / (1 + c_qj * (N + sp.sqrt(N)) / alpha))}

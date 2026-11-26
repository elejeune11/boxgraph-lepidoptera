def run(inputs, settings):
    M, N, kappa = inputs["M"], inputs["N"], inputs["kappa"]
    return {"kappa_r_M": M * (N - M) * kappa / N}


def symbolic():
    import sympy as sp
    M, N, kappa = sp.symbols("M N kappa", positive=True)
    return {"kappa_r_M": M * (N - M) * kappa / N}

import math


def run(inputs, settings):
    R_n = float(inputs["R_n"])
    kon_D = float(inputs["kon_over_D"])
    koff_D = float(inputs["koff_over_D"])
    N_b = float(inputs["N_b"])
    P = R_n * math.exp(2.0 * kon_D * math.sqrt(N_b)) * math.exp(-koff_D * N_b)
    return {"P": P}


def symbolic():
    import sympy as sp

    R_n, kon_over_D, koff_over_D, N_b = sp.symbols(
        "R_n kon_over_D koff_over_D N_b", positive=True
    )
    return {
        "P": R_n * sp.exp(2 * kon_over_D * sp.sqrt(N_b)) * sp.exp(-koff_over_D * N_b)
    }

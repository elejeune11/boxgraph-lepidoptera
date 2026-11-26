def run(inputs, settings):
    return {"L_p": 2.0 * inputs["kappa"] / inputs["T"]}


def symbolic():
    import sympy as sp
    kappa, T = sp.symbols("kappa T", positive=True)
    return {"L_p": 2 * kappa / T}

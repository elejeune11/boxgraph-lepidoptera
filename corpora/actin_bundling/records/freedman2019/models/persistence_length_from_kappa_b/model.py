KB_PN_UM_PER_K = 1.380649e-5


def run(inputs, settings):
    kBT = KB_PN_UM_PER_K * inputs["T"]
    return {"L_p": inputs["kappa_B"] / kBT}


def symbolic():
    import sympy as sp

    kappa_B, T = sp.symbols("kappa_B T", positive=True)
    return {"L_p": kappa_B / (KB_PN_UM_PER_K * T)}

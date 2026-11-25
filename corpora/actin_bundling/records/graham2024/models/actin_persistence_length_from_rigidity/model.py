KB_PN_UM_PER_K = 1.380649e-5  # 1.380649e-23 J/K * 1e18 (pN*um / J)


def run(inputs, settings):
    kBT = KB_PN_UM_PER_K * inputs["T"]
    return {"L_p": inputs["K_bend"] / kBT}


def symbolic():
    import sympy as sp

    K_bend, T = sp.symbols("K_bend T", positive=True)
    return {"L_p": K_bend / (KB_PN_UM_PER_K * T)}

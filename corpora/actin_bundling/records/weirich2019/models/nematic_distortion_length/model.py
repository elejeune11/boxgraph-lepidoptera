def run(inputs, settings):
    # 1 nN/m = 1e-3 pN/um -> K [pN] / (gamma [nN/m] * 1e-3 [pN/um per nN/m]) = um
    gamma_pN_per_um = inputs["gamma"] * 1e-3
    return {"l_K": inputs["K"] / gamma_pN_per_um}


def symbolic():
    import sympy as sp
    K, gamma = sp.symbols("K gamma", positive=True)
    return {"l_K": K / (gamma * sp.Rational(1, 1000))}

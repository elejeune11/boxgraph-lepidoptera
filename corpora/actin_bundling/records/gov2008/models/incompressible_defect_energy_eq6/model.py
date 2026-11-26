def run(inputs, settings):
    K = inputs["K"]
    rho = inputs["rho"]
    R = inputs["R"]
    L_p = inputs["L_p"]
    a = inputs["a"]
    return {"E_defects_b": K * rho * ((R / L_p) ** 4) / ((R / a) ** 2)}


def symbolic():
    import sympy as sp

    K, rho, R, L_p, a = sp.symbols("K rho R L_p a", positive=True)
    return {"E_defects_b": K * rho * (R / L_p) ** 4 / (R / a) ** 2}

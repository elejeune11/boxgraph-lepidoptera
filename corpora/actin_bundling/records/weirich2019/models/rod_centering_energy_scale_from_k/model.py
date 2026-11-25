def run(inputs, settings):
    k = inputs["k"]      # pN/um
    L_D = inputs["L_D"]  # um
    L_M = inputs["L_M"]  # um
    # k [pN/um] * um^2 = pN*um -> * 1000 = pN*nm
    F_0_pN_um = k * ((L_M ** 2 - 4.0 * L_D ** 2) ** 2) / (16.0 * (L_M ** 2 + 4.0 * L_D ** 2))
    return {"F_0": F_0_pN_um * 1000.0}


def symbolic():
    import sympy as sp
    k, L_D, L_M = sp.symbols("k L_D L_M", positive=True)
    return {"F_0": 1000 * k * (L_M ** 2 - 4 * L_D ** 2) ** 2 / (16 * (L_M ** 2 + 4 * L_D ** 2))}

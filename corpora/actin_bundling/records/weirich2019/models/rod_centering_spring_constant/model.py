def run(inputs, settings):
    F_0_pN_um = inputs["F_0"] * 1e-3  # pN*nm -> pN*um
    L_D = inputs["L_D"]               # um
    L_M = inputs["L_M"]               # um
    k_over_F0 = 16.0 * (L_M ** 2 + 4.0 * L_D ** 2) / ((L_M ** 2 - 4.0 * L_D ** 2) ** 2)
    return {
        "k": F_0_pN_um * k_over_F0,
        "k_over_F0": k_over_F0,
    }


def symbolic():
    import sympy as sp
    F_0, L_D, L_M = sp.symbols("F_0 L_D L_M", positive=True)
    k_over_F0 = 16 * (L_M ** 2 + 4 * L_D ** 2) / (L_M ** 2 - 4 * L_D ** 2) ** 2
    return {
        "k": (F_0 * sp.Rational(1, 1000)) * k_over_F0,
        "k_over_F0": k_over_F0,
    }

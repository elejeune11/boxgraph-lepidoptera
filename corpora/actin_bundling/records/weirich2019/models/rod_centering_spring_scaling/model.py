def run(inputs, settings):
    F_0_pN_um = inputs["F_0"] * 1e-3  # pN*nm -> pN*um
    L_D = inputs["L_D"]               # um
    return {"k": F_0_pN_um / (L_D ** 2)}


def symbolic():
    import sympy as sp
    F_0, L_D = sp.symbols("F_0 L_D", positive=True)
    return {"k": (F_0 * sp.Rational(1, 1000)) / L_D ** 2}

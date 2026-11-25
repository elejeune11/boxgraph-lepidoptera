def run(inputs, settings):
    # kBT in pN*nm -> pN*um is * 1e-3
    kBT_pN_um = inputs["kBT"] * 1e-3
    return {"mu": inputs["D"] / kBT_pN_um}


def symbolic():
    import sympy as sp

    D, kBT = sp.symbols("D kBT", positive=True)
    return {"mu": D / (kBT * sp.Rational(1, 1000))}

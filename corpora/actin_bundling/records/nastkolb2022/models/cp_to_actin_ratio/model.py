def run(inputs, settings):
    # c_cp in nM -> uM is * 1e-3; c_actin in uM
    return {"R_cp": (inputs["c_cp"] * 1e-3) / inputs["c_actin"]}


def symbolic():
    import sympy as sp

    c_cp, c_actin = sp.symbols("c_cp c_actin", positive=True)
    return {"R_cp": sp.Rational(1, 1000) * c_cp / c_actin}

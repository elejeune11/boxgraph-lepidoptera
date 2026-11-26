def run(inputs, settings):
    return {"kon_over_D": float(inputs["two_kon_over_D"]) / 2.0}


def symbolic():
    import sympy as sp

    two_kon_over_D = sp.symbols("two_kon_over_D", positive=True)
    return {"kon_over_D": two_kon_over_D / 2}

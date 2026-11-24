def run(inputs, settings):
    return {"E_bend_0": 0.5 * inputs["k_b"] * 1e-3 * inputs["L"] / (inputs["R"] ** 2)}


def symbolic():
    import sympy as sp

    k_b, L, R = sp.symbols("k_b L R", positive=True)
    return {"E_bend_0": sp.Rational(1, 2000) * k_b * L / R ** 2}

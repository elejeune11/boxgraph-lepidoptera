def run(inputs, settings):
    return {"R_mean": (inputs["t"] * inputs["D_B"] * inputs["w"]) ** (1.0 / 3.0)}


def symbolic():
    import sympy as sp

    t, D_B, w = sp.symbols("t D_B w", positive=True)
    return {"R_mean": (t * D_B * w) ** sp.Rational(1, 3)}

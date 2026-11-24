def run(inputs, settings):
    # Literally as printed: RHS (t * D_B / w^2)^(1/3) is dimensionless rather than [length]
    return {"R_mean": (inputs["t"] * inputs["D_B"] / (inputs["w"] ** 2)) ** (1.0 / 3.0)}


def symbolic():
    import sympy as sp

    t, D_B, w = sp.symbols("t D_B w", positive=True)
    return {"R_mean": (t * D_B / w**2) ** sp.Rational(1, 3)}

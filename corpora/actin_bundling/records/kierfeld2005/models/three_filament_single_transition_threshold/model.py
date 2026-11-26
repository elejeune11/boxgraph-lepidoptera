def run(inputs, settings):
    # From 8 * kappa_123 = kappa_12^r with r2 = kappa_2 / kappa_1 and r3 = kappa_3 / kappa_1:
    # 8 * (1 + r2) * r3 / (1 + r2 + r3) = r2 / (1 + r2)
    # => r3_c = r2 * (1 + r2) / (8 * (1 + r2)^2 - r2).
    # For kappa_1 = kappa_2 (r2 = 1): r3_c = 2 / 31 = 0.064516...
    r2 = inputs["r2"]
    return {"r3_c": r2 * (1.0 + r2) / (8.0 * (1.0 + r2) ** 2 - r2)}


def symbolic():
    import sympy as sp
    r2 = sp.symbols("r2", positive=True)
    return {"r3_c": r2 * (1 + r2) / (8 * (1 + r2) ** 2 - r2)}

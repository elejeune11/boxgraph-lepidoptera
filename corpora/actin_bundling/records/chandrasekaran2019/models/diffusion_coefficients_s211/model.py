def run(inputs, settings):
    return {
        "D_alpha": inputs["r_alpha"] * inputs["D_actin"],
        "D_myosin": inputs["r_myosin"] * inputs["D_actin"],
    }


def symbolic():
    import sympy as sp
    D_actin, r_alpha, r_myosin = sp.symbols("D_actin r_alpha r_myosin", positive=True)
    return {
        "D_alpha": r_alpha * D_actin,
        "D_myosin": r_myosin * D_actin,
    }

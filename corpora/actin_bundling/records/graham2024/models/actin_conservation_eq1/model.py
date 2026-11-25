def run(inputs, settings):
    T = inputs["T_actin"]
    F = inputs["F_actin"]
    return {
        "G_actin": T - F,
        "free_fraction": 1.0 - F / T,
    }


def symbolic():
    import sympy as sp

    T_actin, F_actin = sp.symbols("T_actin F_actin", positive=True)
    return {
        "G_actin": T_actin - F_actin,
        "free_fraction": 1 - F_actin / T_actin,
    }

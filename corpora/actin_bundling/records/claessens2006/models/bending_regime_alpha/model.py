def run(inputs, settings):
    return {
        "alpha": (inputs["k_par"] * inputs["L_B"] ** 2)
        / (inputs["delta"] * inputs["k_a"])
    }


def symbolic():
    import sympy as sp

    k_par, L_B, delta, k_a = sp.symbols("k_par L_B delta k_a", positive=True)
    return {"alpha": k_par * L_B**2 / (delta * k_a)}

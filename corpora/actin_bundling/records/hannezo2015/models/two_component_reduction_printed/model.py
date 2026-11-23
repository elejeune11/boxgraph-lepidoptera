def run(inputs, settings):
    K = inputs["K"]
    return {
        "D": inputs["D_m"] / (1.0 + K),
        "k_d_eff": K * inputs["k_d"] / (1.0 + K),
    }


def symbolic():
    import sympy as sp

    D_m, k_d, K = sp.symbols("D_m k_d K", positive=True)
    return {
        "D": D_m / (1 + K),
        "k_d_eff": K * k_d / (1 + K),
    }

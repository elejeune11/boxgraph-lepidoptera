def run(inputs, settings):
    K = inputs["K"]
    return {
        "D": K * inputs["D_m"] / (1.0 + K),
        "k_d_eff": inputs["k_d"] / (1.0 + K),
    }


def symbolic():
    import sympy as sp

    D_m, k_d, K = sp.symbols("D_m k_d K", positive=True)
    return {
        "D": K * D_m / (1 + K),
        "k_d_eff": k_d / (1 + K),
    }

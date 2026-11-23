def run(inputs, settings):
    return {
        "K2_bar": 2.0 * inputs["K2"],
        "K3_bar": 0.25 * inputs["K3"],
    }


def symbolic():
    import sympy as sp
    K2, K3 = sp.symbols("K2 K3", positive=True)
    return {
        "K2_bar": 2 * K2,
        "K3_bar": K3 / 4,
    }

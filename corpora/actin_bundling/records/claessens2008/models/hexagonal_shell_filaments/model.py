def run(inputs, settings):
    k = inputs["k"]
    return {
        "n_f": 1.0 + 3.0 * k * (k + 1.0),
        "D_fil": 2.0 * k + 1.0,
    }


def symbolic():
    import sympy as sp

    k = sp.Symbol("k", positive=True)
    return {
        "n_f": 1 + 3 * k * (k + 1),
        "D_fil": 2 * k + 1,
    }

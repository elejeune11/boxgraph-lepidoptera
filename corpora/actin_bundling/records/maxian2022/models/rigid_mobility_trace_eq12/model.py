def run(inputs, settings):
    a = inputs["alpha"]
    b = inputs["beta"]
    return {
        "c_tr": 3.0 * a + b,
        "anisotropy_ratio": (a + b) / a,
    }


def symbolic():
    import sympy as sp

    alpha, beta = sp.symbols("alpha beta", positive=True)
    return {
        "c_tr": 3 * alpha + beta,
        "anisotropy_ratio": (alpha + beta) / alpha,
    }

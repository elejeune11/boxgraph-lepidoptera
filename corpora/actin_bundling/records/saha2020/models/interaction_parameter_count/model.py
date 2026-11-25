def run(inputs, settings):
    N = float(inputs["N"])
    n_pair = N * (N - 1.0) / 2.0
    return {
        "n_alpha": n_pair,
        "n_chi": n_pair,
        "n_total": N * (N - 1.0),
    }


def symbolic():
    import sympy as sp

    N = sp.Symbol("N", positive=True)
    return {
        "n_alpha": N * (N - 1) / 2,
        "n_chi": N * (N - 1) / 2,
        "n_total": N * (N - 1),
    }

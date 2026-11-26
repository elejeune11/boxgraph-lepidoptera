def run(inputs, settings):
    return {"pct": 100.0 * inputs["n_mode"] / inputs["n_total"]}


def symbolic():
    import sympy as sp

    n_mode, n_total = sp.symbols("n_mode n_total", positive=True)
    return {"pct": 100 * n_mode / n_total}

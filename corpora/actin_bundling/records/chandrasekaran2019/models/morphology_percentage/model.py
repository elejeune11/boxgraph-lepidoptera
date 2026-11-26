def run(inputs, settings):
    return {"pct": 100.0 * inputs["n_class"] / inputs["N_total"]}


def symbolic():
    import sympy as sp
    n_class, N_total = sp.symbols("n_class N_total", positive=True)
    return {"pct": 100 * n_class / N_total}

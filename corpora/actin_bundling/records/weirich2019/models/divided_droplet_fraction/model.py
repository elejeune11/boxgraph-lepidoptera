def run(inputs, settings):
    return {"frac_divided": 100.0 * inputs["N_divided"] / inputs["N_total"]}


def symbolic():
    import sympy as sp
    N_divided, N_total = sp.symbols("N_divided N_total", positive=True)
    return {"frac_divided": 100 * N_divided / N_total}

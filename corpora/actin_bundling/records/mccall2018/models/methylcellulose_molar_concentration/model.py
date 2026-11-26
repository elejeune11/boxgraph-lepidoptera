def run(inputs, settings):
    # c_wv in g/dL (1 g/dL = 10 g/L); M_W in kDa (1 kDa = 1000 g/mol)
    # c [mol/L] = (10 * c_wv) / (1000 * M_W) = 0.01 * c_wv / M_W -> c [uM] = 1e4 * c_wv / M_W
    return {"c": 1.0e4 * inputs["c_wv"] / inputs["M_W"]}


def symbolic():
    import sympy as sp

    c_wv, M_W = sp.symbols("c_wv M_W", positive=True)
    return {"c": 10000 * c_wv / M_W}

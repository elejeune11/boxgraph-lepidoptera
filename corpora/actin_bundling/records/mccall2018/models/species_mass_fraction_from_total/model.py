def run(inputs, settings):
    return {"f_sp": inputs["f_tot"] / 2.0}


def symbolic():
    import sympy as sp

    f_tot = sp.Symbol("f_tot", positive=True)
    return {"f_sp": f_tot / 2}

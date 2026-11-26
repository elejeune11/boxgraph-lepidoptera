def run(inputs, settings):
    f_mem = float(inputs["f_mem"])  # pN
    f_fil = float(inputs["f_fil"])  # pN
    return {"N": f_mem / f_fil}


def symbolic():
    import sympy as sp

    f_mem, f_fil = sp.symbols("f_mem f_fil", positive=True)
    return {"N": f_mem / f_fil}

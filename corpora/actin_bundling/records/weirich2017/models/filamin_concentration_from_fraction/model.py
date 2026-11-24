def run(inputs, settings):
    return {"c_fil": inputs["c_A"] * inputs["f_fil"]}


def symbolic():
    import sympy as sp

    c_A, f_fil = sp.symbols("c_A f_fil", positive=True)
    return {"c_fil": c_A * f_fil}

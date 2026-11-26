def run(inputs, settings):
    kon_D = float(inputs["kon_over_D"])
    koff_D = float(inputs["koff_over_D"])
    return {"N_peak": (kon_D / koff_D) ** 2}


def symbolic():
    import sympy as sp

    kon_over_D, koff_over_D = sp.symbols("kon_over_D koff_over_D", positive=True)
    return {"N_peak": (kon_over_D / koff_over_D) ** 2}

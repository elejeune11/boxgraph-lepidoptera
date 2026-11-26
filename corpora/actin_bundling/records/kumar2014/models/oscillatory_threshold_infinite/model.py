def run(inputs, settings):
    return {"P_c": (1.0 + inputs["alpha"]) / (inputs["A0_fA"] + inputs["I0_fI"])}


def symbolic():
    import sympy as sp

    alpha = sp.symbols("alpha", positive=True)
    A0_fA, I0_fI = sp.symbols("A0_fA I0_fI", real=True)
    return {"P_c": (1 + alpha) / (A0_fA + I0_fI)}

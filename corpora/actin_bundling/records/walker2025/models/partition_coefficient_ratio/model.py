def run(inputs, settings):
    return {"P": float(inputs["I_in"]) / float(inputs["I_out"])}


def symbolic():
    import sympy as sp

    I_in, I_out = sp.symbols("I_in I_out", positive=True)
    return {"P": I_in / I_out}

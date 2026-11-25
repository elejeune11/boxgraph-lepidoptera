def run(inputs, settings):
    return {"P_I": inputs["I_inside"] / inputs["I_outside"]}


def symbolic():
    import sympy as sp

    I_inside, I_outside = sp.symbols("I_inside I_outside", positive=True)
    return {"P_I": I_inside / I_outside}

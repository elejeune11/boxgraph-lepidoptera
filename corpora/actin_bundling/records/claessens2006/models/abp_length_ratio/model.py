def run(inputs, settings):
    return {"ratio": inputs["t_aact"] / inputs["t_fascin"]}


def symbolic():
    import sympy as sp

    t_aact, t_fascin = sp.symbols("t_aact t_fascin", positive=True)
    return {"ratio": t_aact / t_fascin}

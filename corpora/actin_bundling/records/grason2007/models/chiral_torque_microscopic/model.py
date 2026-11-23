def run(inputs, settings):
    return {"gamma": inputs["dE"] / (inputs["l"] * inputs["P"])}


def symbolic():
    import sympy as sp
    dE, l, P = sp.symbols("dE l P", positive=True)
    return {"gamma": dE / (l * P)}

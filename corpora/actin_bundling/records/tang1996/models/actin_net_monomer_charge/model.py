def run(inputs, settings):
    return {"q_net": inputs["q_raw"] - inputs["n_His"]}


def symbolic():
    import sympy as sp

    q_raw, n_His = sp.symbols("q_raw n_His", positive=True)
    return {"q_net": q_raw - n_His}

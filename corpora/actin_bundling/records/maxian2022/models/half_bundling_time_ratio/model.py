def run(inputs, settings):
    return {"ratio": inputs["t50_nb"] / inputs["t50_b"]}


def symbolic():
    import sympy as sp

    t50_nb, t50_b = sp.symbols("t50_nb t50_b", positive=True)
    return {"ratio": t50_nb / t50_b}

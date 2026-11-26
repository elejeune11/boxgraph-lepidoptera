def run(inputs, settings):
    return {"ratio": inputs["t0_nb"] / inputs["t0_b"]}


def symbolic():
    import sympy as sp

    t0_nb, t0_b = sp.symbols("t0_nb t0_b", positive=True)
    return {"ratio": t0_nb / t0_b}

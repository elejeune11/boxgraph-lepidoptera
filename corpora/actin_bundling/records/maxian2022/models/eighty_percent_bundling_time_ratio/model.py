def run(inputs, settings):
    return {"ratio": inputs["t80_long"] / inputs["t80_short"]}


def symbolic():
    import sympy as sp

    t80_long, t80_short = sp.symbols("t80_long t80_short", positive=True)
    return {"ratio": t80_long / t80_short}

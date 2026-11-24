def run(inputs, settings):
    return {"L_ratio": 1.0 + inputs["n_cycles"] * inputs["step_frac"]}


def symbolic():
    import sympy as sp

    n_cycles, step_frac = sp.symbols("n_cycles step_frac", positive=True)
    return {"L_ratio": 1 + n_cycles * step_frac}

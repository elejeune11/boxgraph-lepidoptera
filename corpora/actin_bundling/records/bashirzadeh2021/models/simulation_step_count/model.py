def run(inputs, settings):
    return {"N_steps": inputs["t_end"] / inputs["dt"]}


def symbolic():
    import sympy as sp

    t_end, dt = sp.symbols("t_end dt", positive=True)
    return {"N_steps": t_end / dt}

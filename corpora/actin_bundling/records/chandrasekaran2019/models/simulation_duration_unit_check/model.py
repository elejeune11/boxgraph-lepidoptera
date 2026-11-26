def run(inputs, settings):
    return {"t_min": inputs["t_sec"] / 60.0}


def symbolic():
    import sympy as sp
    t_sec = sp.symbols("t_sec", positive=True)
    return {"t_min": t_sec / 60}

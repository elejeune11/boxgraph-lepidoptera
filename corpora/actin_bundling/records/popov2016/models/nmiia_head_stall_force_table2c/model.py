def run(inputs, settings):
    return {"F_stall": inputs["K_head"] * inputs["d_step"]}

def symbolic():
    import sympy as sp
    K_head, d_step = sp.symbols("K_head d_step", positive=True)
    return {"F_stall": K_head * d_step}

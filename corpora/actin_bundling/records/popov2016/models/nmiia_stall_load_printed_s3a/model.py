def run(inputs, settings):
    return {"F_ext": inputs["rho"] * inputs["K_head"] * inputs["d_step"]}

def symbolic():
    import sympy as sp
    rho, K_head, d_step = sp.symbols("rho K_head d_step", positive=True)
    return {"F_ext": rho * K_head * d_step}

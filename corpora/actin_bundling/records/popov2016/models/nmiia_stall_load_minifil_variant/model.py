def run(inputs, settings):
    return {"F_ext": inputs["rho"] * inputs["N_total"] * inputs["K_head"] * inputs["d_step"]}

def symbolic():
    import sympy as sp
    rho, N_total, K_head, d_step = sp.symbols("rho N_total K_head d_step", positive=True)
    return {"F_ext": rho * N_total * K_head * d_step}

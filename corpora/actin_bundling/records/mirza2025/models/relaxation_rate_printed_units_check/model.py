def run(inputs, settings):
    return {"ratio": inputs["K_S_printed"] / inputs["eta_rot"]}

def symbolic():
    import sympy as sp
    K_S_printed, eta_rot = sp.symbols("K_S_printed eta_rot", positive=True)
    return {"ratio": K_S_printed / eta_rot}

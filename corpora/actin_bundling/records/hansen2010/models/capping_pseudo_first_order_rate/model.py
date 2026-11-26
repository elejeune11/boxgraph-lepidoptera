def run(inputs, settings):
    return {"K_obs": inputs["K_plus"] * inputs["C_cp"]}


def symbolic():
    import sympy as sp

    K_plus, C_cp = sp.symbols("K_plus C_cp", positive=True)
    return {"K_obs": K_plus * C_cp}

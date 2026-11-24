def run(inputs, settings):
    # L_target, L0 in um; k_grow in nm/s -> 1e-3 um/s
    return {"t_grow": (inputs["L_target"] - inputs["L0"]) / (inputs["k_grow"] * 1e-3)}


def symbolic():
    import sympy as sp

    L_target, L0, k_grow = sp.symbols("L_target L0 k_grow", positive=True)
    return {"t_grow": 1000 * (L_target - L0) / k_grow}

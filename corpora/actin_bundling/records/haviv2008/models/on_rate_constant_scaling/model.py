def run(inputs, settings):
    # Constant on-rate proportional to bulk actin density n_a (independent of phi)
    _ = float(inputs["phi"])
    return {"two_kon_over_D": float(inputs["n_a"])}


def symbolic():
    import sympy as sp

    n_a, phi = sp.symbols("n_a phi", positive=True)
    return {"two_kon_over_D": n_a}

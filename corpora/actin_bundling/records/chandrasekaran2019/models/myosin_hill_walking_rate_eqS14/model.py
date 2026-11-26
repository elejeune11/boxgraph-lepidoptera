def run(inputs, settings):
    fs = inputs["F_stall"]
    fe = inputs["F_ext"]
    z = inputs["zeta"]
    return {"k_walk": inputs["k_walk_0"] * (fs - fe) / (fs + fe / z)}


def symbolic():
    import sympy as sp
    k_walk_0, F_stall, F_ext, zeta = sp.symbols("k_walk_0 F_stall F_ext zeta", positive=True)
    return {"k_walk": k_walk_0 * (F_stall - F_ext) / (F_stall + F_ext / zeta)}

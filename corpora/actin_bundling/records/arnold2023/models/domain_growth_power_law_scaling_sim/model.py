"""Power-law coarsening scaling a ~ t^alpha in dimensionless simulation time (prefactor=unknown)."""


def run(inputs, settings):
    t = float(inputs["t"])
    alpha = float(inputs["alpha"])
    val = t ** alpha
    return {"a_norm": val, "a_sim": val}


def symbolic():
    import sympy as sp

    t, alpha = sp.symbols("t alpha", positive=True)
    return {"a_norm": t ** alpha, "a_sim": t ** alpha}

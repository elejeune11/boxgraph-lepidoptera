"""Power-law coarsening scaling a/a_0 ~ t^alpha (prefactor=unknown; prefactor set to 1)."""


def run(inputs, settings):
    t = float(inputs["t"])
    alpha = float(inputs["alpha"])
    return {"a_norm": t ** alpha}


def symbolic():
    import sympy as sp

    t, alpha = sp.symbols("t alpha", positive=True)
    return {"a_norm": t ** alpha}

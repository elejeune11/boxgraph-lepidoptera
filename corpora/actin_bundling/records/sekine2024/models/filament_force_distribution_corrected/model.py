"""Corrected force assignment matching Fig. S1 and force/torque conservation (PDF p. 29)."""


def run(inputs, settings):
    f = inputs["f"]
    alpha = inputs["alpha"]
    f_0 = (1.0 - alpha) * f
    f_1 = alpha * f
    return {"f_0": f_0, "f_1": f_1, "f_sum_ratio": (f_0 + f_1) / f}


def symbolic():
    import sympy as sp

    f, alpha = sp.symbols("f alpha", positive=True)
    f_0 = (1 - alpha) * f
    f_1 = alpha * f
    return {"f_0": f_0, "f_1": f_1, "f_sum_ratio": sp.Integer(1)}

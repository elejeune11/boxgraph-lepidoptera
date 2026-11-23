"""Printed force assignment to filament endpoint particles (SI Notes p. 2, PDF p. 29).

Printed form has a minus sign before delta_{i1}:
  f_0 = (1 - alpha) * f
  f_1 = -alpha * f
so (f_0 + f_1) / f = 1 - 2 * alpha != 1 for alpha != 0.
"""


def run(inputs, settings):
    f = inputs["f"]
    alpha = inputs["alpha"]
    f_0 = (1.0 - alpha) * f
    f_1 = -alpha * f
    return {"f_0": f_0, "f_1": f_1, "f_sum_ratio": (f_0 + f_1) / f}


def symbolic():
    import sympy as sp

    f, alpha = sp.symbols("f alpha", positive=True)
    f_0 = (1 - alpha) * f
    f_1 = -alpha * f
    return {"f_0": f_0, "f_1": f_1, "f_sum_ratio": 1 - 2 * alpha}

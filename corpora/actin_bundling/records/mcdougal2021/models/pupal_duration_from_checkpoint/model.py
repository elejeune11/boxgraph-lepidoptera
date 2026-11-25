"""Implied total pupal duration T_pupal = 100 * t / p_dev (inverse of p. 6, SI p. 10)."""


def run(inputs, settings):
    return {"T_pupal": 100.0 * float(inputs["t"]) / float(inputs["p_dev"])}


def symbolic():
    import sympy as sp

    t, p_dev = sp.symbols("t p_dev", positive=True)
    return {"T_pupal": 100 * t / p_dev}

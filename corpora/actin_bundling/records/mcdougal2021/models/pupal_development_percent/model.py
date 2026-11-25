"""Normalized pupal development percentage p_dev = 100 * t / T_pupal (p. 6, SI p. 10)."""


def run(inputs, settings):
    return {"p_dev": 100.0 * float(inputs["t"]) / float(inputs["T_pupal"])}


def symbolic():
    import sympy as sp

    t, T_pupal = sp.symbols("t T_pupal", positive=True)
    return {"p_dev": 100 * t / T_pupal}

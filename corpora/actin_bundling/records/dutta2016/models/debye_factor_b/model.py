import math


def run(inputs, settings):
    k = abs(float(inputs["k"]))
    if k < 1e-2:
        k2 = k * k
        # Exact Taylor expansion around k=0 to avoid 0/0 catastrophic cancellation
        val = 1.0 - k2 / 6.0 + (k2 * k2) / 48.0 - (k2 * k2 * k2) / 480.0
    else:
        val = (-8.0 + 8.0 * math.exp(-0.5 * k * k) + 4.0 * k * k) / (k ** 4)
    return {"b": val}


def symbolic():
    import sympy as sp
    k = sp.Symbol("k", positive=True)
    return {"b": (-8 + 8 * sp.exp(-k**2 / 2) + 4 * k**2) / k**4}

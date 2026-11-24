import math


def _binom_half_sq(n):
    val = 1.0
    for k in range(n):
        val *= (0.5 - k) / (k + 1)
    return val * val


def run(inputs, settings):
    a = inputs["a"]
    b = inputs["b"]
    n_terms = max(1, int(round(settings.get("n_terms", 7))))
    h = ((a - b) / (a + b)) ** 2
    s = 0.0
    hp = 1.0
    for n in range(n_terms):
        s += _binom_half_sq(n) * hp
        hp *= h
    return {"L": math.pi * (a + b) * s}


def symbolic():
    import sympy as sp

    a, b = sp.symbols("a b", positive=True)
    h = ((a - b) / (a + b)) ** 2
    s = sum(sp.binomial(sp.Rational(1, 2), n) ** 2 * h ** n for n in range(7))
    return {"L": sp.pi * (a + b) * s}

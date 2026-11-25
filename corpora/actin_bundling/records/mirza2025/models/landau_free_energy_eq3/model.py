def run(inputs, settings):
    # L [pN = Pa*um^2] * grad_q^2 [1/um^2] = Pa
    S = inputs["S"]
    gq = inputs["grad_q"]
    f = 0.5 * inputs["a"] * S**2 + 0.125 * inputs["b"] * S**4 + 0.5 * inputs["L"] * gq**2
    return {"f": f}

def symbolic():
    import sympy as sp
    a, b, L, S, grad_q = sp.symbols("a b L S grad_q", real=True)
    return {"f": sp.Rational(1, 2) * a * S**2 + sp.Rational(1, 8) * b * S**4 + sp.Rational(1, 2) * L * grad_q**2}

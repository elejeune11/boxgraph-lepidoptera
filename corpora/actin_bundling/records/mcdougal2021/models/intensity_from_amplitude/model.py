"""Reflected intensity I = A^2 (SI Fig. S3B caption, p. 15)."""


def run(inputs, settings):
    A = float(inputs["A"])
    return {"I": A * A}


def symbolic():
    import sympy as sp

    A = sp.Symbol("A", real=True)
    return {"I": A**2}

"""Approximate rectangular wing scale area A = L * W (SI p. 5)."""


def run(inputs, settings):
    return {"A": float(inputs["L"]) * float(inputs["W"])}


def symbolic():
    import sympy as sp

    L, W = sp.symbols("L W", positive=True)
    return {"A": L * W}

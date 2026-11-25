import math


def run(inputs, settings):
    decay = math.exp(-inputs["t"] / inputs["tau"])
    A = inputs["A_f"] + (inputs["A_i"] - inputs["A_f"]) * decay
    return {"A": A, "A_norm": decay}


def symbolic():
    import sympy as sp

    A_i, A_f, t, tau = sp.symbols("A_i A_f t tau", positive=True)
    decay = sp.exp(-t / tau)
    return {"A": A_f + (A_i - A_f) * decay, "A_norm": decay}

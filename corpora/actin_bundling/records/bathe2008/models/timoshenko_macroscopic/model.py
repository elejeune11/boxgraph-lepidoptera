"""Macroscopic Timoshenko beam three-point bending stiffness (p. 5 col 1)."""


def run(inputs, settings):
    E_B = inputs["E_B"]
    I_B = inputs["I_B"]
    beta = inputs["beta"]
    G_B = inputs["G_B"]
    A_B = inputs["A_B"]
    L = inputs["L"]
    return {"kappa_B": (E_B * I_B) / (1.0 + 12.0 * E_B * I_B / (beta * G_B * A_B * (L ** 2)))}


def symbolic():
    import sympy as sp

    E_B, I_B, beta, G_B, A_B, L = sp.symbols("E_B I_B beta G_B A_B L", positive=True)
    return {"kappa_B": (E_B * I_B) / (1 + 12 * E_B * I_B / (beta * G_B * A_B * L**2))}

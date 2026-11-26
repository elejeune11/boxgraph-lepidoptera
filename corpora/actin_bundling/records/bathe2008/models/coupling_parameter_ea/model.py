"""Fiber coupling parameter alpha = (k_x / (E_f * A_f)) * (L^2 / delta)."""


def run(inputs, settings):
    return {"alpha": (inputs["k_x"] / inputs["E_f_A_f"]) * (inputs["L"] ** 2 / inputs["delta"])}


def symbolic():
    import sympy as sp

    k_x, E_f_A_f, L, delta = sp.symbols("k_x E_f_A_f L delta", positive=True)
    return {"alpha": (k_x / E_f_A_f) * (L**2 / delta)}

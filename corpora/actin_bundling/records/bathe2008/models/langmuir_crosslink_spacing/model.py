"""Langmuir isotherm mean cross-link spacing (Fig. 4 caption): delta = delta_min * (K_d + c_fascin) / c_fascin."""


def run(inputs, settings):
    return {"delta": inputs["delta_min"] * (inputs["K_d"] + inputs["c_fascin"]) / inputs["c_fascin"]}


def symbolic():
    import sympy as sp

    delta_min, K_d, c_fascin = sp.symbols("delta_min K_d c_fascin", positive=True)
    return {"delta": delta_min * (K_d + c_fascin) / c_fascin}

def run(inputs, settings):
    return {"k_b_c": 16.0 * inputs["k_0"] * (inputs["phi_0_B"] ** 3)}


def symbolic():
    import sympy as sp

    k_0, phi_0_B = sp.symbols("k_0 phi_0_B", positive=True)
    return {"k_b_c": 16 * k_0 * phi_0_B**3}

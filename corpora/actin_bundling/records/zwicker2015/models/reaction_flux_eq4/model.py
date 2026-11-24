def run(inputs, settings):
    return {
        "s": inputs["k_f"] * inputs["phi_A"]
        - inputs["k_b"] * inputs["phi_B"]
        + inputs["k"] * inputs["phi_A"] * inputs["phi_B"]
    }


def symbolic():
    import sympy as sp

    k_f, k_b, k, phi_A, phi_B = sp.symbols("k_f k_b k phi_A phi_B", positive=True)
    return {"s": k_f * phi_A - k_b * phi_B + k * phi_A * phi_B}

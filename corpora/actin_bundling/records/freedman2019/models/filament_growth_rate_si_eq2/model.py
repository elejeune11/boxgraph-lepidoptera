def run(inputs, settings):
    return {
        "k_grow": inputs["k_g"] * inputs["delta_l_a"] * 1000.0,
        "P_step": inputs["k_g"] * inputs["dt"],
    }


def symbolic():
    import sympy as sp

    k_g, delta_l_a, dt = sp.symbols("k_g delta_l_a dt", positive=True)
    return {
        "k_grow": 1000 * k_g * delta_l_a,
        "P_step": k_g * dt,
    }

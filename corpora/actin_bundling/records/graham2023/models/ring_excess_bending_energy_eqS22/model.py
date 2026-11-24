def run(inputs, settings):
    return {
        "Delta_E_bend": 0.5 * inputs["k_b"] * 1e-3 * (inputs["L"] / (inputs["R"] ** 4))
        * (inputs["u_m"] ** 2) * (2.0 + inputs["m"] ** 4)
    }


def symbolic():
    import sympy as sp

    k_b, L, R, m, u_m = sp.symbols("k_b L R m u_m", positive=True)
    return {"Delta_E_bend": sp.Rational(1, 2000) * k_b * (L / R ** 4) * u_m ** 2 * (2 + m ** 4)}

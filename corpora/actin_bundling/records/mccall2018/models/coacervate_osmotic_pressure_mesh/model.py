K_B = 1.380649e-23  # J/K (exact SI)


def run(inputs, settings):
    xi_m = inputs["xi"] * 1.0e-9
    return {"Pi": (K_B * inputs["T"]) / (xi_m ** 3)}


def symbolic():
    import sympy as sp

    xi, T = sp.symbols("xi T", positive=True)
    return {"Pi": sp.Float(K_B * 1.0e27) * T / (xi ** 3)}

import math


def run(inputs, settings):
    q_c = inputs["q_c"]
    # First non-uniform Neumann (no-flux) mode on [0, L]: q_1 = pi / L <= q_c
    L_c = math.pi / q_c
    lambda_c = 2.0 * math.pi / q_c
    return {
        "L_c": L_c,
        "lambda_c": lambda_c,
        "wavelength_to_length_ratio": lambda_c / L_c,
    }


def symbolic():
    import sympy as sp

    q_c = sp.Symbol("q_c", positive=True)
    L_c = sp.pi / q_c
    lambda_c = 2 * sp.pi / q_c
    return {
        "L_c": L_c,
        "lambda_c": lambda_c,
        "wavelength_to_length_ratio": sp.simplify(lambda_c / L_c),
    }

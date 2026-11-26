import math

E_CHARGE = 1.602176634e-19       # C (exact SI)
EPSILON_0 = 8.8541878128e-12     # F/m (CODATA)
K_B = 1.380649e-23               # J/K (exact SI)


def run(inputs, settings):
    eps = inputs["epsilon"]
    T = inputs["T"]
    lb_m = (E_CHARGE ** 2) / (4.0 * math.pi * EPSILON_0 * eps * K_B * T)
    return {"lambda_B": lb_m * 1.0e10}


def symbolic():
    import sympy as sp

    epsilon, T = sp.symbols("epsilon T", positive=True)
    coeff = (E_CHARGE ** 2) * 1.0e10 / (4.0 * math.pi * EPSILON_0 * K_B)
    return {"lambda_B": sp.Float(coeff) / (epsilon * T)}

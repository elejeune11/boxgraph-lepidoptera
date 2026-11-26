import math


def run(inputs, settings):
    # Printed Eq. (11): dP/dt = (1/2)*D_R*d^2P/dR^2 + d(R_dot*P)/dR = d/dR[(1/2)*D_R*P' + R_dot*P].
    # With R_dot = -beta * dE_total/dR (Eq. 10), zero flux (1/2)*D_R*P' - beta*(dE/dR)*P = 0
    # integrates to P(R) ~ exp[+2 * beta * (E_total - E_0) / D_R], i.e. boltzmann_coeff = +2.
    drift_sign = +1.0  # printed '+' in Eq. (11)
    boltzmann_coeff = 2.0 * drift_sign
    beta = inputs["beta"]
    E_total = inputs["E_total"]
    E_0 = inputs["E_0"]
    D_R = inputs["D_R"]
    return {
        "boltzmann_coeff": boltzmann_coeff,
        "P": math.exp(boltzmann_coeff * beta * (E_total - E_0) / D_R),
    }


def symbolic():
    import sympy as sp

    beta, E_total, E_0, D_R = sp.symbols("beta E_total E_0 D_R", positive=True)
    return {
        "boltzmann_coeff": sp.Integer(2),
        "P": sp.exp(2 * beta * (E_total - E_0) / D_R),
    }

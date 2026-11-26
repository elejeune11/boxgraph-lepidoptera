import math


def run(inputs, settings):
    eps_0 = float(inputs["eps_0"])
    U = float(inputs["U"])
    C = float(inputs["C"])
    mu = float(inputs["mu"])
    kBT = float(inputs["kBT"])
    w0 = float(inputs["omega_0"])
    beta = 1.0 / kBT
    xi_t = 2.0 * beta * C
    c_U = 0.5
    eps_0_prime = eps_0 - c_U * U
    z = math.exp(beta * (mu + eps_0_prime))
    # Printed sites: (0, 0) -> (4, 1) has dl = 4, dm = +1; (4, 1) -> (9, -1) has dl = 5, dm = -2
    V_nn_4 = 0.5 * beta * (U ** 2) * math.cos(w0 * 4.0 - 2.0 * math.pi * 1.0 / 6.0) * math.exp(-4.0 / xi_t)
    V_nn_5 = 0.5 * beta * (U ** 2) * math.cos(w0 * 5.0 - 2.0 * math.pi * (-2.0) / 6.0) * math.exp(-5.0 / xi_t)
    return {
        "xi_t": xi_t,
        "eps_0_prime": eps_0_prime,
        "c_U": c_U,
        "z": z,
        "V_nn_4": V_nn_4,
        "V_nn_5": V_nn_5,
    }


def symbolic():
    import sympy as sp

    eps_0, U, C, mu, kBT, omega_0 = sp.symbols("eps_0 U C mu kBT omega_0", positive=True)
    beta = 1 / kBT
    xi_t = 2 * beta * C
    eps_0_prime = eps_0 - U / 2
    return {
        "xi_t": xi_t,
        "eps_0_prime": eps_0_prime,
        "c_U": sp.Rational(1, 2),
        "z": sp.exp(beta * (mu + eps_0_prime)),
        "V_nn_4": sp.Rational(1, 2) * beta * U ** 2 * sp.cos(4 * omega_0 - sp.pi / 3) * sp.exp(-4 / xi_t),
        "V_nn_5": sp.Rational(1, 2) * beta * U ** 2 * sp.cos(5 * omega_0 + 2 * sp.pi / 3) * sp.exp(-5 / xi_t),
    }

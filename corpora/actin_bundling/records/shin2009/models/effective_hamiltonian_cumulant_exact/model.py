import math


def run(inputs, settings):
    eps_0 = float(inputs["eps_0"])
    U = float(inputs["U"])
    C = float(inputs["C"])
    kBT = float(inputs["kBT"])
    w0 = float(inputs["omega_0"])
    beta = 1.0 / kBT
    # Under H_0 = (C/2)*sum(Delta_phi_l - w0)^2, Delta_phi_l ~ N(w0, 1/(beta*C)) i.i.d.
    # => <exp(i*(phi_l - phi_l'))>_0 = exp(i*w0*(l - l')) * exp(-|l - l'|/(2*beta*C))
    # => xi_t = 2*beta*C
    xi_t = 2.0 * beta * C
    # First-order cumulant of -eps_0 + U*(1 - cos(phi_l - phi_m)):
    # <cos(phi_l - phi_m)>_0 = 0 (uniform zero mode) => -eps_0 + U = -(eps_0 - U), so c_U = 1.0
    c_U_exact = 1.0
    eps_0_prime_O1 = eps_0 - c_U_exact * U
    # Second-order cumulant:
    # (1/2)*(beta*U)^2 * <cos(phi_l - phi_m)*cos(phi_l' - phi_m')>_0
    # = (1/2)*beta * [(beta*U^2/2) * cos(w0*(l-l') - (phi_m - phi_m')) * exp(-|l-l'|/xi_t)]
    # Self term (l = l'): n_l^2 = n_l contributes +(1/2)*(beta*U^2/2) = +beta*U^2/4 to eps_0'
    eps_0_prime_O2_no_self = eps_0_prime_O1 + 0.25 * beta * (U ** 2)
    V_nn_4 = 0.5 * beta * (U ** 2) * math.cos(w0 * 4.0 - 2.0 * math.pi * 1.0 / 6.0) * math.exp(-4.0 / xi_t)
    V_nn_5 = 0.5 * beta * (U ** 2) * math.cos(w0 * 5.0 - 2.0 * math.pi * (-2.0) / 6.0) * math.exp(-5.0 / xi_t)
    return {
        "xi_t": xi_t,
        "eps_0_prime_O1": eps_0_prime_O1,
        "eps_0_prime_O2_no_self": eps_0_prime_O2_no_self,
        "c_U_exact": c_U_exact,
        "V_nn_4": V_nn_4,
        "V_nn_5": V_nn_5,
    }


def symbolic():
    import sympy as sp

    eps_0, U, C, kBT, omega_0 = sp.symbols("eps_0 U C kBT omega_0", positive=True)
    beta = 1 / kBT
    xi_t = 2 * beta * C
    return {
        "xi_t": xi_t,
        "eps_0_prime_O1": eps_0 - U,
        "eps_0_prime_O2_no_self": eps_0 - U + sp.Rational(1, 4) * beta * U ** 2,
        "c_U_exact": sp.Integer(1),
        "V_nn_4": sp.Rational(1, 2) * beta * U ** 2 * sp.cos(4 * omega_0 - sp.pi / 3) * sp.exp(-4 / xi_t),
        "V_nn_5": sp.Rational(1, 2) * beta * U ** 2 * sp.cos(5 * omega_0 + 2 * sp.pi / 3) * sp.exp(-5 / xi_t),
    }

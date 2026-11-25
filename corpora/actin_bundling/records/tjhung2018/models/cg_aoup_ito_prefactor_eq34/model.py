def run(inputs, settings):
    # For Ito SDE with noise sqrt(2*D(r)) where D = tau*v^2 + T (grad D = 2*tau*v*grad v)
    # and drift a = s * tau*v*grad v = (s/2) * grad D:
    # J = (a - grad D)*rho - D*grad rho = -rho*D * grad[log rho + (1 - s/2)*log D].
    # Printed Eq. (34) has s = -1 (-tau*v*grad_i v), yielding prefactor 1 - (-1)/2 = 3/2.
    s = -1.0
    return {"log_prefactor": 1.0 - 0.5 * s}


def symbolic():
    import sympy as sp

    tau, v0 = sp.symbols("tau v0", positive=True)
    return {"log_prefactor": sp.Rational(3, 2) + 0 * tau * v0}

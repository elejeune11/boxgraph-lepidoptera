def run(inputs, settings):
    c = float(inputs["c"])
    alpha = float(inputs["alpha"])
    kappa = float(inputs["kappa"])
    q = float(inputs["q"])
    q2 = q * q
    # Printed Eq. (9): partial_t A = nabla^2[-(c^2 + i*alpha/2)*A + ...] + kappa*nabla^4 A
    # For A ~ exp(lambda*t + i*q*x): nabla^2 -> -q^2, nabla^4 -> +q^4
    # => lambda = (c^2 + i*alpha/2)*q^2 + kappa*q^4 (UV-unstable +kappa*q^4 and rescaled coefficients!)
    return {
        "re_lambda": c * c * q2 + kappa * q2 * q2,
        "im_lambda": 0.5 * alpha * q2,
    }


def symbolic():
    import sympy as sp

    c, alpha, kappa, q = sp.symbols("c alpha kappa q", real=True)
    return {
        "re_lambda": c**2 * q**2 + kappa * q**4,
        "im_lambda": alpha * q**2 / 2,
    }

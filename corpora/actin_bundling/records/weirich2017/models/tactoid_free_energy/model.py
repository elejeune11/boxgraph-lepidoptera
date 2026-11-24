import math


def run(inputs, settings):
    R = inputs["R"]
    a = inputs["alpha"]
    f = math.sin(a) - a * math.cos(a)
    g = f - (1.0 / 3.0) * (math.sin(a) ** 3)
    A_um2 = 4.0 * math.pi * (R**2) * f
    A_omega_um2 = 2.0 * math.pi * (R**2) * g
    # (nN/m) * (um^2) = (1e-9 N/m) * (1e-12 m^2) = 1e-21 J
    F_J = (inputs["gamma_I"] * A_um2 + inputs["gamma_A"] * A_omega_um2) * 1e-21
    return {"F": F_J}


def symbolic():
    import sympy as sp

    R, alpha, gamma_I, gamma_A = sp.symbols("R alpha gamma_I gamma_A", positive=True)
    f = sp.sin(alpha) - alpha * sp.cos(alpha)
    g = f - sp.Rational(1, 3) * sp.sin(alpha) ** 3
    A = 4 * sp.pi * R**2 * f
    A_omega = 2 * sp.pi * R**2 * g
    return {"F": (gamma_I * A + gamma_A * A_omega) * sp.Rational(1, 10**21)}

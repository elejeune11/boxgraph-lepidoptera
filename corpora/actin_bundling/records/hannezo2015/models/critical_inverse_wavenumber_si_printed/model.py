import math


def run(inputs, settings):
    # Printed on page-011.png: 1/q_c = pi * L * (D * tau * xi / eta)^(1/4)
    # Note: D*tau*xi/eta is dimensionless: (m^2/s * s * Pa*s/m) / (Pa*s*m) = 1
    D_si = inputs["D"] * 1e-12
    ratio = D_si * inputs["tau"] * inputs["xi"] / inputs["eta"]
    return {"inv_qc": math.pi * inputs["L"] * (ratio ** 0.25)}


def symbolic():
    import sympy as sp

    L, D, tau, xi, eta = sp.symbols("L D tau xi eta", positive=True)
    D_si = D * sp.Rational(1, 10**12)
    return {"inv_qc": sp.pi * L * (D_si * tau * xi / eta) ** sp.Rational(1, 4)}

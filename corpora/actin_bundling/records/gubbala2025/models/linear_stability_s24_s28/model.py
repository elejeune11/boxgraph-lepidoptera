import math


def run(inputs, settings):
    gamma = float(inputs["gamma"])
    chi = float(inputs["chi"])
    beta = float(inputs["beta"])
    E1 = float(inputs["E1"])
    E3 = float(inputs["E3"])
    kx = float(inputs["kx"])
    ky = float(inputs["ky"])

    k2 = kx * kx + ky * ky
    sb = math.sqrt(beta)

    # Printed Eqs. (S27), (S28), (S26)
    omega_phi = gamma * k2 * (2.0 + (1.0 + 1.5 * chi) * kx * kx + (1.0 - 0.5 * chi) * ky * ky)
    omega_xx = 4.0 + (E1 + E3 / sb) * kx * kx + (E1 - E3 / sb) * ky * ky
    omega_xy = E1 * k2

    # Printed 2x2 Jacobian from Eqs. (S24)-(S25):
    # d(phi~)/dt = -omega_phi * phi~ + (4*gamma/sqrt(beta))*k^2 * Qxx~
    # d(Qxx~)/dt = (1/sqrt(beta)) * phi~ - omega_xx * Qxx~
    c12 = (4.0 * gamma / sb) * k2
    c21 = 1.0 / sb
    tr = -(omega_phi + omega_xx)
    det = omega_phi * omega_xx - c12 * c21
    disc = tr * tr - 4.0 * det
    sqrt_disc = math.sqrt(max(0.0, disc))
    sigma_max = 0.5 * (tr + sqrt_disc)
    sigma_min = 0.5 * (tr - sqrt_disc)

    # Printed high-k positivity bounds from Eqs. (S27)-(S28)
    chi_min = -2.0 / 3.0
    chi_max = 2.0
    beta_min = (E3 / E1) ** 2

    return {
        "omega_phi": omega_phi,
        "omega_xx": omega_xx,
        "omega_xy": omega_xy,
        "sigma_max": sigma_max,
        "sigma_min": sigma_min,
        "chi_min": chi_min,
        "chi_max": chi_max,
        "beta_min": beta_min,
    }


def symbolic():
    import sympy as sp

    gamma, beta, E1 = sp.symbols("gamma beta E1", positive=True)
    chi, E3, kx, ky = sp.symbols("chi E3 kx ky", real=True)
    k2 = kx**2 + ky**2
    sb = sp.sqrt(beta)
    omega_phi = gamma * k2 * (2 + (1 + sp.Rational(3, 2) * chi) * kx**2 + (1 - chi / 2) * ky**2)
    omega_xx = 4 + (E1 + E3 / sb) * kx**2 + (E1 - E3 / sb) * ky**2
    omega_xy = E1 * k2
    c12 = 4 * gamma * k2 / sb
    c21 = 1 / sb
    tr = -(omega_phi + omega_xx)
    det = omega_phi * omega_xx - c12 * c21
    sqrt_disc = sp.sqrt(tr**2 - 4 * det)
    return {
        "omega_phi": omega_phi,
        "omega_xx": omega_xx,
        "omega_xy": omega_xy,
        "sigma_max": (tr + sqrt_disc) / 2,
        "sigma_min": (tr - sqrt_disc) / 2,
        "chi_min": -sp.Rational(2, 3),
        "chi_max": sp.Integer(2),
        "beta_min": (E3 / E1) ** 2,
    }

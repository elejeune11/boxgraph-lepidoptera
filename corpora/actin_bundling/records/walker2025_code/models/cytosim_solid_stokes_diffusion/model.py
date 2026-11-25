"""Port of CytoSim Solid Stokes drag diffusivity (solid.cc:854-943) with Wrist diffusion no-op."""

import math


def run(inputs, settings):
    eta = float(inputs["eta"])
    r = float(inputs["r_solid"])
    kT = float(inputs["kT"])
    dt = float(inputs["dt"])

    drag = 6.0 * math.pi * eta * r  # pN * s / um
    D_solid = kT / drag  # um^2 / s
    sigma_B_nm = math.sqrt(2.0 * D_solid * dt) * 1000.0  # nm
    return {"D_solid": D_solid, "sigma_B": sigma_B_nm}


def symbolic():
    import sympy as sp

    eta, r_solid, kT, dt = sp.symbols("eta r_solid kT dt", positive=True)
    D_solid = kT / (6 * sp.pi * eta * r_solid)
    sigma_B = sp.sqrt(2 * D_solid * dt) * 1000
    return {"D_solid": D_solid, "sigma_B": sigma_B}

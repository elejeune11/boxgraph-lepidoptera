"""Faithful Python port of CytoSim Sphere/Solid/Bead Stokes drag vs Couple uniform diffusion:
- src/sim/sphere.cc:275-284 (spDrag = 6*pi*eta*R, spDragRot = 8*pi*eta*R^3)
- src/sim/solid.cc:1289-1292 (soDrag = 6*pi*sumR, soDragRot = 8*pi*sumR3; scaled by |viscosity|)
- src/sim/bead.cc:128-136 (paDrag = 6*pi*eta*R unless bead:drag > 0)
- src/sim/couple_prop.cc:146, src/sim/single_prop.cc:207 (diffusion_dt = sqrt(6*D*dt))
"""
import math


def run(inputs, settings):
    r_um = float(inputs["radius"]) * 1e-3
    eta = abs(float(inputs["viscosity"]))
    kT = float(inputs["kT"])
    dt = float(inputs["dt"])
    d_presc = float(inputs["D_prescribed"])

    vr = eta * r_um
    drag_t = (6.0 * math.pi) * vr
    drag_r = (8.0 * math.pi) * (vr * r_um * r_um)
    mob_stokes = 1.0 / drag_t
    d_stokes = kT * mob_stokes
    sigma_stokes = math.sqrt(2.0 * d_stokes * dt)

    hw_couple = math.sqrt(6.0 * d_presc * dt)  # couple_prop.cc:146
    sigma_couple = hw_couple / math.sqrt(3.0)  # == sqrt(2 * d_presc * dt)
    eta_req = kT / (6.0 * math.pi * r_um * d_presc)

    return {
        "drag_trans": drag_t,
        "drag_rot": drag_r,
        "mobility_stokes": mob_stokes,
        "D_stokes": d_stokes,
        "sigma_gauss_stokes": sigma_stokes,
        "halfwidth_uniform_couple": hw_couple,
        "sigma_uniform_couple": sigma_couple,
        "eta_required_for_D": eta_req,
    }


def symbolic():
    import sympy as sp

    radius, viscosity, kT, dt, D_prescribed = sp.symbols(
        "radius viscosity kT dt D_prescribed", positive=True
    )
    r_um = radius * sp.Rational(1, 1000)
    drag_t = 6 * sp.pi * viscosity * r_um
    drag_r = 8 * sp.pi * viscosity * r_um ** 3
    mob = 1 / drag_t
    d_stokes = kT / drag_t
    return {
        "drag_trans": drag_t,
        "drag_rot": drag_r,
        "mobility_stokes": mob,
        "D_stokes": d_stokes,
        "sigma_gauss_stokes": sp.sqrt(2 * d_stokes * dt),
        "halfwidth_uniform_couple": sp.sqrt(6 * D_prescribed * dt),
        "sigma_uniform_couple": sp.sqrt(2 * D_prescribed * dt),
        "eta_required_for_D": kT / (6 * sp.pi * r_um * D_prescribed),
    }

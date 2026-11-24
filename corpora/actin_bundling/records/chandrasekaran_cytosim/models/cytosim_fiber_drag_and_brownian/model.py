"""Faithful Python port of CytoSim cylindrical fiber drag & Brownian scaling:
- src/sim/chain.cc:36-40 (Chain::bestNumberOfPoints)
- src/sim/fiber.cc:884, 976-1002 (Fiber::dragCoefficientCylinder)
- src/sim/fiber.cc:1057-1084 (Fiber::setDragCoefficient)
- src/sim/mecafil.cc:89-98 (Mecafil::addBrownianForces)
- src/sim/meca.cc:346, 662 (alpha = 2 * kT / dt; noiseLevel *= tau)
"""
import math


def _best_number_of_points(ratio: float) -> int:
    n = 1.0 + math.floor(ratio)
    return int(n + (1 if ratio * (n - 0.5) >= n * (n - 1.0) else 0))


def run(inputs, settings):
    L = float(inputs["L"])
    seg = float(inputs["segmentation"])
    visc = abs(float(inputs["viscosity"]))
    visc_water = float(inputs["viscosity_water"])
    r = float(inputs["drag_radius"])
    l_cut = float(inputs["drag_length"])
    kT = float(inputs["kT"])
    dt = float(inputs["dt"])

    nbp = _best_number_of_points(L / seg)

    # Port of Fiber::dragCoefficientCylinder (src/sim/fiber.cc:976-1002)
    DRAG = 3.0  # src/sim/fiber.cc:884
    drag_coeff = 6.0 * r
    lenc = max(min(L, l_cut), r)
    drag_cylinder = DRAG * L / (math.log(0.5 * lenc / r) + 0.32)
    min_len = 2.0 * r * math.exp(1.0 - 0.32)
    if L > min_len:
        drag_coeff = max(drag_coeff, drag_cylinder)
    drag_total = math.pi * visc * drag_coeff

    # Port of Fiber::setDragCoefficient (src/sim/fiber.cc:1077)
    point_mob = nbp / drag_total
    com_mob = 1.0 / drag_total
    d_com = kT * com_mob

    # Port of Mecafil::addBrownianForces (src/sim/mecafil.cc:91-97) + Meca::solve (meca.cc:346,662)
    alpha = 2.0 * kT / dt
    b = math.sqrt(alpha / point_mob)
    disp_vertex = dt * b * point_mob  # == sqrt(2 * kT * dt * point_mob)
    disp_com = math.sqrt(2.0 * d_com * dt)

    return {
        "N_pts": float(nbp),
        "drag_total": drag_total,
        "point_mobility": point_mob,
        "com_mobility": com_mob,
        "D_com": d_com,
        "brownian_disp_vertex": disp_vertex,
        "brownian_disp_com": disp_com,
        "viscosity_ratio_water": visc / visc_water,
    }

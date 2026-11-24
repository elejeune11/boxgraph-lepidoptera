"""Faithful Python port of CytoSim discrete fiber bending mechanics:
- src/sim/chain.cc:36-40 (Chain::bestNumberOfPoints)
- src/sim/chain.cc:1170-1193 (Chain::bendingEnergy0)
- src/sim/fiber.h:257 (Fiber::bendingEnergy)
- src/sim/fiber.cc:1097 (Fiber::prepareMecable: iRigidity = rigidity / s^3)
- src/sim/mecafil.cc:189-203, 276-313 (add_rigidity0 / add_rigidity3D)
- src/sim/fiber_prop.h:59 (rigidity = Lp * kT)
"""
import math


def _best_number_of_points(ratio: float) -> int:
    # Port of Chain::bestNumberOfPoints (src/sim/chain.cc:36-40)
    n = 1.0 + math.floor(ratio)
    return int(n + (1 if ratio * (n - 0.5) >= n * (n - 1.0) else 0))


def _add_rigidity_3d(pts, r1: float):
    # Port of add_rigidity0 / add_rigidity3D (src/sim/mecafil.cc:189-203, 276-313)
    nbp = len(pts)
    forces = [[0.0, 0.0, 0.0] for _ in range(nbp)]
    r2 = 2.0 * r1
    for p in range(nbp - 2):
        for d in range(3):
            f = (pts[p][d] + pts[p + 2][d]) - 2.0 * pts[p + 1][d]
            forces[p][d] -= f * r1
            forces[p + 1][d] += f * r2
            forces[p + 2][d] -= f * r1
    return forces


def run(inputs, settings):
    rigidity = float(inputs["rigidity"])
    L = float(inputs["L"])
    seg_target = float(inputs["segmentation"])
    R = float(inputs["R_curve"])
    kT = float(inputs["kT"])

    nbp = _best_number_of_points(L / seg_target)
    n_seg = nbp - 1
    s = L / n_seg
    i_rigidity = rigidity / (s ** 3)

    # Construct an exact discrete circular arc of chord length s on a circle of radius R
    # (sin(dtheta / 2) = s / (2 * R))
    dtheta = 2.0 * math.asin(min(1.0, s / (2.0 * R)))
    pts = [
        [R * math.cos(k * dtheta), R * math.sin(k * dtheta), 0.0]
        for k in range(nbp)
    ]

    # Exact quadratic Hamiltonian differentiated by add_rigidity3D
    e_ham = 0.0
    lsp = nbp - 2
    dot_sum = 0.0
    for p in range(lsp):
        ax, ay, az = pts[p]
        bx, by, bz = pts[p + 1]
        cx, cy, cz = pts[p + 2]
        fx = ax + cx - 2.0 * bx
        fy = ay + cy - 2.0 * by
        fz = az + cz - 2.0 * bz
        e_ham += 0.5 * i_rigidity * (fx * fx + fy * fy + fz * fz)
        dot_sum += (bx - ax) * (cx - bx) + (by - ay) * (cy - by) + (bz - az) * (cz - bz)

    # Port of Chain::bendingEnergy0 (src/sim/chain.cc:1170-1193) * rigidity (src/sim/fiber.h:257)
    if lsp > 0:
        e0 = lsp - dot_sum / (s * s)
        e0 *= (lsp + 1.0) / (s * lsp)
        e_cyto = rigidity * e0
    else:
        e_cyto = 0.0

    forces = _add_rigidity_3d(pts, i_rigidity)
    max_f = max(math.sqrt(f[0] ** 2 + f[1] ** 2 + f[2] ** 2) for f in forces)
    e_cont = 0.5 * rigidity * L / (R * R)
    lp = rigidity / kT

    return {
        "N_pts": float(nbp),
        "N_seg": float(n_seg),
        "s_actual": s,
        "iRigidity": i_rigidity,
        "E_hamiltonian": e_ham,
        "E_cyto_reported": e_cyto,
        "E_continuum": e_cont,
        "max_node_force": max_f,
        "persistence_length": lp,
    }

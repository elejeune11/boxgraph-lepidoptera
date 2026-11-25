"""Derived quantities implied directly by walker2025_code .cym configs and CytoSim default kT."""

import math

AVOGADRO = 6.02214076e23


def run(inputs, settings):
    n_steps = float(inputs["n_steps"])
    nb_frames = float(inputs["nb_frames"])
    dt = float(inputs["dt"])
    L0 = float(inputs["L0"])
    k_grow = float(inputs["k_grow"])
    L_max = float(inputs["L_max"])
    N_blob = float(inputs["N_blob"])
    r_solid = float(inputs["r_solid"])
    R_drop = float(inputs["R_drop"])
    rigidity = float(inputs["rigidity"])
    kT = float(inputs["kT"])

    T_sim = n_steps * dt
    dt_frame = T_sim / nb_frames
    n_snaps_out = nb_frames + 1.0
    L_final = min(L0 + k_grow * T_sim, L_max)

    # Sphere volume in liters: 1 um^3 = 1e-15 L; mol/L -> uM multiplies by 1e6
    V_L = (4.0 / 3.0) * math.pi * (R_drop ** 3) * 1e-15
    C_uM = (N_blob / AVOGADRO) / V_L * 1e6

    phi_pct = N_blob * ((r_solid / R_drop) ** 3) * 100.0
    L_p = rigidity / kT

    return {
        "T_sim": T_sim,
        "dt_frame": dt_frame,
        "n_snaps_out": n_snaps_out,
        "L_final": L_final,
        "C_crosslinker": C_uM,
        "phi_crosslinker": phi_pct,
        "L_p": L_p,
    }

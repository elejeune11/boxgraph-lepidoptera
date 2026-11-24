"""Faithful Python port of CytoSim monomer-depletion and Brownian ratchet polymerization:
- src/sim/fiber_set.cc:226-249 (FiberSet::steps: free_polymer = max(0, 1 - used_polymer / total_polymer))
- src/sim/fibers/growing_fiber.cc:66-79 (GrowingFiber::step: one-sided exp(forceP / growing_force) for forceP < 0)
- src/sim/fibers/classic_fiber.cc:112-125 (ClassicFiber::stepPlusEnd: identical growth rule)
- src/sim/chain.cc:1673-1676 (Chain::projectedForceEndP: negative when antagonizing elongation)
"""
import math


def run(inputs, settings):
    k_grow_nm_s = float(inputs["growing_speed"])
    f_grow = float(inputs["growing_force"])
    f_proj = float(inputs["projected_force"])
    used_poly = float(inputs["used_polymer"])
    total_poly = float(inputs["total_polymer"])
    n_fil = float(inputs["n_filaments"])
    L0_um = float(inputs["L0"])
    R_um = float(inputs["R_drop"])

    # Port of FiberSet::steps() (src/sim/fiber_set.cc:239-248)
    free_poly = 1.0 - used_poly / total_poly
    if free_poly < 0.0:
        free_poly = 0.0

    # Port of GrowingFiber::step() (src/sim/fibers/growing_fiber.cc:72-76)
    # Note: one-sided clamp! Only antagonistic force (f_proj < 0) multiplies by exp(f_proj / f_grow).
    v_fil = k_grow_nm_s * free_poly
    if f_proj < 0.0 and v_fil > 0.0:
        v_fil *= math.exp(f_proj / f_grow)

    v_net = n_fil * v_fil
    v_fil_um_s = v_fil * 1e-3
    t_to_2r = (2.0 * R_um - L0_um) / v_fil_um_s if v_fil_um_s > 0.0 else float("inf")

    return {
        "free_polymer": free_poly,
        "v_per_filament": v_fil,
        "v_total_network": v_net,
        "t_to_2R": t_to_2r,
    }

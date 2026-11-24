"""Faithful Python port of CytoSim coupled Nucleator + FiberSet monomer depletion:
- src/sim/hands/nucleator_prop.cc:33 (default hold_end = MINUS_END)
- src/sim/hands/nucleator.cc:140-144, 171-201 (Nucleator::stepUnattached: creates fiber & attaches minus end)
- src/sim/single.cc:148-151 (Single::stepA: detaches at unbinding_rate to become free again)
- src/sim/fiber_set.cc:215-249 (FiberSet::steps: free_polymer = max(0, 1 - used_polymer / total_polymer))
- src/sim/fibers/growing_fiber.cc:66-79 (GrowingFiber::step: addP = growing_speed * dt * free_polymer)
- src/sim/fiber.cc:180-203 (Fiber::updateLength: caps each filament at max_length)

Zero antagonistic end force (force-free bulk growth limit) is evaluated; exact mean-field
step evolution matches the ensemble average of CytoSim's Gillespie nucleator cycle.
"""
import math

AVOGADRO = 6.02214076e23


def run(inputs, settings):
    n_seed = float(inputs["N_seed"])
    l_seed_init = float(inputs["L_seed_init"])
    l_max = float(inputs["L_max"])
    l_nuc_init = float(inputs["L_nuc_init"])
    k_grow_um_s = float(inputs["k_grow"]) * 1e-3
    c_nuc_nM = float(inputs["C_nucleator"])
    k_nuc = float(inputs["k_nuc"])
    k_off = float(inputs["k_unbind_nuc"])
    c_total_uM = float(inputs["C_total_actin"])
    r_drop_um = float(inputs["R_drop"])
    dl_um = float(inputs["dl_monomer"]) * 1e-3
    t_eval = float(inputs["t_eval"])

    dt = float(settings.get("dt", 0.02))

    # Droplet volume in liters: 1 um^3 = 1e-15 L
    vol_L = (4.0 * math.pi / 3.0) * (r_drop_um ** 3) * 1e-15
    n_nuc = float(round(c_nuc_nM * 1e-9 * AVOGADRO * vol_L))
    total_polymer = (c_total_uM * 1e-6 * AVOGADRO * vol_L) * dl_um

    # Exact discrete-step probabilities matching CytoSim's Simul::steps() order:
    # 1. singles.steps() -> Nucleator::stepUnattached / Single::stepA (simul_step.cc:85)
    # 2. fibers.steps()  -> FiberSet::steps() updates free_polymer, then GrowingFiber::step() (simul_step.cc:90)
    p_nuc = -math.expm1(-k_nuc * dt)
    p_off = -math.expm1(-k_off * dt)

    n_steps = max(1, int(round(t_eval / dt)))
    dt_eff = t_eval / n_steps
    if abs(dt_eff - dt) > 1e-12:
        p_nuc = -math.expm1(-k_nuc * dt_eff)
        p_off = -math.expm1(-k_off * dt_eff)

    n_free = n_nuc
    n_bound = 0.0
    n_fil = n_seed
    l_seed = l_seed_init
    used = n_seed * l_seed_init
    # Track total length of nucleated cohort (all nucleated fibers have length <= l_seed)
    for _ in range(n_steps):
        # Step 1: singles.steps() — nucleators fire & bound nucleators detach
        dn_nuc = n_free * p_nuc
        dn_off = n_bound * p_off
        n_free = n_free - dn_nuc + dn_off
        n_bound = n_bound + dn_nuc - dn_off
        n_fil += dn_nuc
        used += dn_nuc * l_nuc_init

        # Step 2: fibers.steps() — update free_polymer (fiber_set.cc:239-248)
        fp = 1.0 - used / total_polymer
        if fp < 0.0:
            fp = 0.0

        # Step 3: GrowingFiber::step() + Fiber::updateLength() (growing_fiber.cc:72, fiber.cc:180-202)
        dl_step = k_grow_um_s * fp * dt_eff
        if l_seed < l_max:
            dl_seed = min(dl_step, l_max - l_seed)
            l_seed += dl_seed
            used += n_seed * dl_seed
        # Newly nucleated filaments are always shorter than seed filaments (< l_max)
        n_nucleated = n_fil - n_seed
        if n_nucleated > 0.0:
            used += n_nucleated * dl_step

    fp_final = max(0.0, 1.0 - used / total_polymer)

    return {
        "N_nuc_molecules": n_nuc,
        "total_polymer_length": total_polymer,
        "N_fil_final": n_fil,
        "L_seed_final": l_seed,
        "free_polymer_final": fp_final,
    }

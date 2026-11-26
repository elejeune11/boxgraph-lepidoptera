"""Faithful Python port of AFINES C++ derived parameters.

Ports:
- prog/bundles.cpp:224, 325-326, 333, 362, 378
- src/motor_ensemble.cpp:42, 313-317 (spacer_ensemble::set_bending: kb = modulus / mld)
- src/bead.cpp:23-25 (friction = 6 * pi * visc * rad)
- src/spacer.cpp:55, 151 (damp = 6 * pi * vis * mld)
- xlink_sorting_paper/README.md:35-37, 201-203 & plots.nb:8466 (L_p = lkb / 0.004)
"""

import math


def run(inputs, settings):
    nmonomer = float(inputs["nmonomer"])
    link_length = float(inputs["link_length"])
    nlink_max = float(inputs["nlink_max"])
    kappa_b = float(inputs["polymer_bending_modulus"])
    kT = float(inputs["temperature"])
    # 1 Pa*s = 1 N*s/m^2 = 1e12 pN * s / 1e12 um^2 = 1 pN*s/um^2
    visc = float(inputs["viscosity"])
    rad_bead = float(inputs["actin_length"])
    l1 = float(inputs["spacer1_length"])
    l2 = float(inputs["spacer2_length"])
    s1_bend = float(inputs["s1_bend"])
    s2_bend = float(inputs["s2_bend"])
    xrange = float(inputs["xrange"])
    yrange = float(inputs["yrange"])
    grid_factor = float(inputs["grid_factor"])
    rho2 = float(inputs["spacer2_density"])
    nfixed = float(inputs["nfixed"])
    mu_lat = float(inputs["lattice_mu"])

    n_links_init = nmonomer - 1.0
    init_len = n_links_init * link_length
    max_len = nlink_max * link_length
    Lp = kappa_b / kT

    # src/motor_ensemble.cpp:314: double kb = modulus / mld;
    kb_short_eff = s1_bend / l1
    kb_long_eff = s2_bend / l2
    kb_ratio_eff = kb_long_eff / kb_short_eff
    s_bend_ratio_raw = s2_bend / s1_bend

    # src/bead.cpp:25 & src/spacer.cpp:55
    damp_bead = 6.0 * math.pi * visc * rad_bead
    mu_bead = 1.0 / damp_bead
    D_bead = kT * mu_bead
    damp_short = 6.0 * math.pi * visc * l1
    D_short = kT / damp_short
    damp_long = 6.0 * math.pi * visc * l2
    D_long = kT / damp_long

    # prog/bundles.cpp:325-326 & src/motor_ensemble.cpp:42
    xgrid = float(int(round(grid_factor * xrange)))
    ygrid = float(int(round(grid_factor * yrange)))
    n_tot = float(int(math.ceil(rho2 * xrange * yrange)))
    n_rand = n_tot - nfixed

    mu_kBT = mu_lat / kT

    return {
        "n_links_init": n_links_init,
        "init_filament_length": init_len,
        "max_filament_length": max_len,
        "persistence_length": Lp,
        "kb_short_eff": kb_short_eff,
        "kb_long_eff": kb_long_eff,
        "kb_ratio_eff": kb_ratio_eff,
        "s_bend_ratio_raw": s_bend_ratio_raw,
        "damp_bead": damp_bead,
        "mu_bead": mu_bead,
        "D_bead": D_bead,
        "damp_short": damp_short,
        "D_short": D_short,
        "damp_long": damp_long,
        "D_long": D_long,
        "xgrid": xgrid,
        "ygrid": ygrid,
        "n_spacers_total": n_tot,
        "n_spacers_random": n_rand,
        "lattice_mu_kBT": mu_kBT,
    }

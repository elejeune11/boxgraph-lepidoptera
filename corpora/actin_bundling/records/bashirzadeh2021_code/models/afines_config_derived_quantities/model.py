"""Faithful port of derived simulation quantities in:
- prog/network.cpp:379-380, 449, 460
- src/bead.cpp:20-22
- src/filament.cpp:67
- src/motor.cpp:35-36
- crosslinker_sorting_in_GUVs/aster_initial_condition_grid_conda_freeze.py:43-45, 54-55, 68-69
"""

import math


def run(inputs, settings):
    xrange = float(inputs["xrange"])
    grid_factor = float(inputs["grid_factor"])
    r_c = float(inputs["r_c"])
    fil_dens = float(inputs["fil_dens"])
    n_beads = float(inputs["n_beads"])
    l_link = float(inputs["l_link"])
    r_bead = float(inputs["r_bead"])
    kappa_b = float(inputs["kappa_b"])
    temperature = float(inputs["temperature"])
    viscosity = float(inputs["viscosity"])  # 1 mg/(um*s) = 1 pN*s/um^2
    dt = float(inputs["dt"])
    tfinal = float(inputs["tfinal"])
    nframes = float(inputs["nframes"])
    rho_fascin = float(inputs["rho_fascin"])
    rho_alpha = float(inputs["rho_alpha"])
    l_fascin = float(inputs["l_fascin"])
    l_alpha = float(inputs["l_alpha"])
    kon_fascin = float(inputs["kon_fascin"])
    koff_fascin = float(inputs["koff_fascin"])
    scl = float(inputs["scl"])
    kon_alpha = float(inputs["kon_alpha"])
    koff_alpha = float(inputs["koff_alpha"])

    area = math.pi * (r_c ** 2)
    n_fil_from_dens = float(round(area * fil_dens))
    n_fascin = float(round(rho_fascin * area))
    n_alpha = float(round(rho_alpha * area))

    L_contour = (n_beads - 1.0) * l_link
    k_b_link = kappa_b / l_link
    L_p_code = kappa_b / temperature

    n_steps = tfinal / dt
    n_bw_print = float(max(int(tfinal / (dt * nframes)), 1))
    dt_frame = tfinal / nframes
    xgrid = float(round(grid_factor * xrange))

    K_fascin = (kon_fascin / scl) / (koff_fascin / scl)
    K_alpha = kon_alpha / koff_alpha
    size_ratio = l_alpha / l_fascin

    dF_bind_fascin = -temperature * math.log(K_fascin)
    dF_bind_alpha = -temperature * math.log(K_alpha)

    # Stokes friction: 6 * pi * vis * r (src/bead.cpp:22, src/motor.cpp:35)
    gamma_bead = 6.0 * math.pi * viscosity * r_bead
    D_bead = temperature / gamma_bead
    gamma_fascin = 6.0 * math.pi * viscosity * l_fascin
    D_fascin = temperature / gamma_fascin
    gamma_alpha = 6.0 * math.pi * viscosity * l_alpha
    D_alpha = temperature / gamma_alpha

    # Leimkuhler-Matthews velocity prefactor: sqrt(T / (2 * damp * dt)) (filament.cpp:67, motor.cpp:36)
    bd_prefactor_bead = math.sqrt(temperature / (2.0 * gamma_bead * dt))
    bd_prefactor_fascin = math.sqrt(temperature / (2.0 * gamma_fascin * dt))
    bd_prefactor_alpha = math.sqrt(temperature / (2.0 * gamma_alpha * dt))

    return {
        "n_fil_from_dens": n_fil_from_dens,
        "n_fascin": n_fascin,
        "n_alpha": n_alpha,
        "L_contour": L_contour,
        "k_b_link": k_b_link,
        "L_p_code": L_p_code,
        "n_steps": n_steps,
        "n_bw_print": n_bw_print,
        "dt_frame": dt_frame,
        "xgrid": xgrid,
        "K_fascin": K_fascin,
        "K_alpha": K_alpha,
        "size_ratio": size_ratio,
        "dF_bind_fascin": dF_bind_fascin,
        "dF_bind_alpha": dF_bind_alpha,
        "gamma_bead": gamma_bead,
        "D_bead": D_bead,
        "gamma_fascin": gamma_fascin,
        "D_fascin": D_fascin,
        "gamma_alpha": gamma_alpha,
        "D_alpha": D_alpha,
        "bd_prefactor_bead": bd_prefactor_bead,
        "bd_prefactor_fascin": bd_prefactor_fascin,
        "bd_prefactor_alpha": bd_prefactor_alpha,
    }

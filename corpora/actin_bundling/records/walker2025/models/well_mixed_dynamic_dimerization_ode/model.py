import math

NA_UM3_PER_UMOL = 602.214076


def run(inputs, settings):
    C_tot = float(inputs["C_tot"])
    C_actin = float(inputs["C_actin"])
    k_bind = float(inputs["k_bind"])
    k_unbind = float(inputs["k_unbind"])
    r_a_um = float(inputs["r_bind_actin"]) * 1e-3
    L_seg = float(inputs["L_seg"])
    k_d_form = float(inputs["k_dimer_form"])
    k_d_split = float(inputs["k_dimer_split"])
    r_m_um = float(inputs["r_bind_monomer"]) * 1e-3

    dt = float(settings.get("dt", 0.002))
    t_end = float(settings.get("t_end", 10.0))

    V_cyl = math.pi * r_a_um * r_a_um * L_seg
    V_sph = (4.0 / 3.0) * math.pi * (r_m_um ** 3)

    # Pseudo-first-order actin binding rate (1/s) and second-order dimerization rate (uM^-1 s^-1)
    k_a_on = k_bind * NA_UM3_PER_UMOL * V_cyl * C_actin
    k_a_off = k_unbind
    k_d_on = k_d_form * NA_UM3_PER_UMOL * V_sph
    k_d_off = k_d_split

    def rhs(state):
        mf, mb, d0, d1, d2 = state
        # Actin binding fluxes
        j_ma = k_a_on * mf - k_a_off * mb
        j_d01 = 2.0 * k_a_on * d0 - k_a_off * d1
        j_d12 = k_a_on * d1 - 2.0 * k_a_off * d2
        # Dimerization fluxes (exact detailed balance across all 3 cycles of Fig. 5A)
        j_ff = 0.5 * k_d_on * mf * mf - k_d_off * d0
        j_fb = k_d_on * mf * mb - k_d_off * d1
        j_bb = 0.5 * k_d_on * mb * mb - k_d_off * d2

        dmf = -j_ma - 2.0 * j_ff - j_fb
        dmb = j_ma - j_fb - 2.0 * j_bb
        dd0 = j_ff - j_d01
        dd1 = j_fb + j_d01 - j_d12
        dd2 = j_bb + j_d12
        return (dmf, dmb, dd0, dd1, dd2)

    n_steps = max(1, int(round(t_end / dt)))
    y = (C_tot, 0.0, 0.0, 0.0, 0.0)
    for _ in range(n_steps):
        k1 = rhs(y)
        y2 = tuple(y[i] + 0.5 * dt * k1[i] for i in range(5))
        k2 = rhs(y2)
        y3 = tuple(y[i] + 0.5 * dt * k2[i] for i in range(5))
        k3 = rhs(y3)
        y4 = tuple(y[i] + dt * k3[i] for i in range(5))
        k4 = rhs(y4)
        y = tuple(
            max(0.0, y[i] + (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]))
            for i in range(5)
        )

    mf, mb, d0, d1, d2 = y
    return {
        "f_free_monomer": mf / C_tot,
        "f_bound_monomer": mb / C_tot,
        "f_free_dimer": 2.0 * d0 / C_tot,
        "f_singly_bound_dimer": 2.0 * d1 / C_tot,
        "f_doubly_bound_dimer": 2.0 * d2 / C_tot,
    }

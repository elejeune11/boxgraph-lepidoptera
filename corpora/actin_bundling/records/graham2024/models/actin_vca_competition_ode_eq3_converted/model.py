import math
from scipy.integrate import solve_ivp

UM_UM3_TO_MOLECULES = 602.214076


def run(inputs, settings):
    V_um3 = (4.0 / 3.0) * math.pi * (inputs["R_drop"] ** 3)
    um_per_uM = V_um3 * UM_UM3_TO_MOLECULES * (inputs["l_mono"] * 1e-3)

    k_grow_uM = inputs["k_grow"] / um_per_uM
    k_bind_tot_uM = (inputs["k_bind"] * inputs["N_VCA"]) / um_per_uM
    N_fil = inputs["N_fil"]
    F_0 = inputs["F_0"]
    T_actin = inputs["T_actin"]
    t_end = inputs["t_end"]
    G_max = inputs["factor_Gmax"] * inputs["C_vca"]
    dt = float(settings.get("dt", 1.0))

    sol0 = solve_ivp(
        lambda t, y: [k_grow_uM * N_fil * (1.0 - y[0] / T_actin)],
        (0.0, t_end),
        [F_0],
        method="RK45",
        max_step=dt,
        rtol=1e-9,
        atol=1e-11,
    )
    F_no_vca = float(sol0.y[0, -1])

    if G_max <= 0.0 or k_bind_tot_uM <= 0.0:
        return {"F_end": F_no_vca, "V_end": 0.0, "delta_F": 0.0}

    def ev_sat(t, y):
        return G_max - y[1]

    ev_sat.terminal = True
    ev_sat.direction = -1

    sol1 = solve_ivp(
        lambda t, y: [
            k_grow_uM * N_fil * (1.0 - (y[0] + y[1]) / T_actin),
            k_bind_tot_uM * (1.0 - (y[0] + y[1]) / T_actin),
        ],
        (0.0, t_end),
        [F_0, 0.0],
        method="RK45",
        events=ev_sat,
        max_step=min(dt, max(1e-6, G_max / max(k_bind_tot_uM, 1.0))),
        rtol=1e-9,
        atol=1e-11,
    )
    t_sat = float(sol1.t[-1])
    F_sat = float(sol1.y[0, -1])
    V_sat = float(sol1.y[1, -1])

    if t_sat >= t_end - 1e-12:
        return {"F_end": F_sat, "V_end": V_sat, "delta_F": F_no_vca - F_sat}

    sol2 = solve_ivp(
        lambda t, y: [k_grow_uM * N_fil * (1.0 - (y[0] + G_max) / T_actin)],
        (t_sat, t_end),
        [F_sat],
        method="RK45",
        max_step=dt,
        rtol=1e-9,
        atol=1e-11,
    )
    F_end = float(sol2.y[0, -1])
    return {"F_end": F_end, "V_end": G_max, "delta_F": F_no_vca - F_end}

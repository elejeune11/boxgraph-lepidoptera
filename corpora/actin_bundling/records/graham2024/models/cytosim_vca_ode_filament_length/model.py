import math
from scipy.integrate import solve_ivp

UM_UM3_TO_MOLECULES = 602.214076


def run(inputs, settings):
    V_um3 = (4.0 / 3.0) * math.pi * (inputs["R_drop"] ** 3)
    N_VCA = math.floor(inputs["C_vca"] * V_um3 * UM_UM3_TO_MOLECULES)
    L_Gmax = (inputs["l_mono"] * 1e-3) * inputs["factor_Gmax"] * N_VCA
    k_grow = inputs["k_grow"]
    N_fil = inputs["N_fil"]
    L_0 = inputs["L_0"]
    L_total = inputs["L_total"]
    t_end = inputs["t_end"]
    k_bind_fiber = inputs["k_bind"] * N_VCA
    dt = float(settings.get("dt", 1.0))

    if L_Gmax <= 0.0 or k_bind_fiber <= 0.0:
        sol = solve_ivp(
            lambda t, y: [k_grow * (1.0 - (N_fil * y[0]) / L_total)],
            (0.0, t_end),
            [L_0],
            method="RK45",
            max_step=dt,
            rtol=1e-9,
            atol=1e-11,
        )
        return {"L_fil": float(sol.y[0, -1])}

    def ev_sat(t, y):
        return L_Gmax - y[1]

    ev_sat.terminal = True
    ev_sat.direction = -1

    def ev_depl(t, y):
        return (1.0 - (N_fil * y[0] + y[1]) / L_total) - 1e-12

    ev_depl.terminal = True
    ev_depl.direction = -1

    jac_mat = [
        [-k_grow * N_fil / L_total, -k_grow / L_total],
        [-k_bind_fiber * N_fil / L_total, -k_bind_fiber / L_total],
    ]

    sol1 = solve_ivp(
        lambda t, y: [
            k_grow * (1.0 - (N_fil * y[0] + y[1]) / L_total),
            k_bind_fiber * (1.0 - (N_fil * y[0] + y[1]) / L_total),
        ],
        (0.0, t_end),
        [L_0, 0.0],
        method="Radau",
        jac=lambda t, y: jac_mat,
        events=[ev_sat, ev_depl],
        first_step=min(dt, max(1e-7, 0.1 * min(L_Gmax, L_total) / max(k_bind_fiber, 1.0))),
        max_step=dt,
        rtol=1e-9,
        atol=1e-11,
    )
    t_sat = float(sol1.t[-1])
    L_sat = float(sol1.y[0, -1])

    if t_sat >= t_end - 1e-12 or len(sol1.t_events[1]) > 0:
        return {"L_fil": L_sat}

    sol2 = solve_ivp(
        lambda t, y: [k_grow * (1.0 - (N_fil * y[0] + L_Gmax) / L_total)],
        (t_sat, t_end),
        [L_sat],
        method="RK45",
        max_step=dt,
        rtol=1e-9,
        atol=1e-11,
    )
    return {"L_fil": float(sol2.y[0, -1])}

from scipy.integrate import solve_ivp


def _integrate_pair(k_grow, N_fil, F_0, T_actin, t_end, k_bind_tot, G_max, dt):
    sol0 = solve_ivp(
        lambda t, y: [k_grow * N_fil * (1.0 - y[0] / T_actin)],
        (0.0, t_end),
        [F_0],
        method="RK45",
        max_step=dt,
        rtol=1e-9,
        atol=1e-11,
    )
    F_no_vca = float(sol0.y[0, -1])

    if G_max <= 0.0 or k_bind_tot <= 0.0:
        return F_no_vca, 0.0, 0.0

    def ev_sat(t, y):
        return G_max - y[1]

    ev_sat.terminal = True
    ev_sat.direction = -1

    sol1 = solve_ivp(
        lambda t, y: [
            k_grow * N_fil * (1.0 - (y[0] + y[1]) / T_actin),
            k_bind_tot * (1.0 - (y[0] + y[1]) / T_actin),
        ],
        (0.0, t_end),
        [F_0, 0.0],
        method="RK45",
        events=ev_sat,
        max_step=min(dt, max(1e-6, G_max / max(k_bind_tot, 1.0))),
        rtol=1e-9,
        atol=1e-11,
    )
    t_sat = float(sol1.t[-1])
    F_sat = float(sol1.y[0, -1])
    V_sat = float(sol1.y[1, -1])

    if t_sat >= t_end - 1e-12:
        return F_sat, V_sat, F_no_vca - F_sat

    sol2 = solve_ivp(
        lambda t, y: [k_grow * N_fil * (1.0 - (y[0] + G_max) / T_actin)],
        (t_sat, t_end),
        [F_sat],
        method="RK45",
        max_step=dt,
        rtol=1e-9,
        atol=1e-11,
    )
    F_end = float(sol2.y[0, -1])
    return F_end, G_max, F_no_vca - F_end


def run(inputs, settings):
    # Evaluates printed Eqs. (1)-(3) literally without converting um/s speeds to uM/s.
    dt = float(settings.get("dt", 1.0))
    G_max = inputs["factor_Gmax"] * inputs["C_vca"]
    k_bind_tot = inputs["k_bind"] * inputs["N_VCA"]
    F_end, V_end, delta_F = _integrate_pair(
        inputs["k_grow"],
        inputs["N_fil"],
        inputs["F_0"],
        inputs["T_actin"],
        inputs["t_end"],
        k_bind_tot,
        G_max,
        dt,
    )
    return {"F_end": F_end, "V_end": V_end, "delta_F": delta_F}

from scipy.integrate import solve_ivp


def run(inputs, settings):
    # Literally evaluates printed SI Eq. (2): d[F-actin]/dt = k_grow * N_fil * (1 - [F-actin]/[T-actin])
    # without converting k_grow (um/s) to concentration units (uM/s).
    k_grow = inputs["k_grow"]
    N_fil = inputs["N_fil"]
    F_0 = inputs["F_0"]
    T_actin = inputs["T_actin"]
    t_end = inputs["t_end"]
    dt = float(settings.get("dt", 1.0))

    def rhs(t, y):
        return [k_grow * N_fil * (1.0 - y[0] / T_actin)]

    sol = solve_ivp(rhs, (0.0, t_end), [F_0], method="RK45", max_step=dt, rtol=1e-9, atol=1e-11)
    F_end = float(sol.y[0, -1])
    return {"F_end": F_end, "G_end": T_actin - F_end}

import math
from scipy.integrate import solve_ivp

UM_UM3_TO_MOLECULES = 602.214076


def run(inputs, settings):
    V_um3 = (4.0 / 3.0) * math.pi * (inputs["R_drop"] ** 3)
    um_per_uM = V_um3 * UM_UM3_TO_MOLECULES * (inputs["l_mono"] * 1e-3)
    k_grow_uM_s = inputs["k_grow"] / um_per_uM

    N_fil = inputs["N_fil"]
    F_0 = inputs["F_0"]
    T_actin = inputs["T_actin"]
    t_end = inputs["t_end"]
    dt = float(settings.get("dt", 1.0))

    def rhs(t, y):
        return [k_grow_uM_s * N_fil * (1.0 - y[0] / T_actin)]

    sol = solve_ivp(rhs, (0.0, t_end), [F_0], method="RK45", max_step=dt, rtol=1e-9, atol=1e-11)
    F_end = float(sol.y[0, -1])
    return {"F_end": F_end, "G_end": T_actin - F_end}

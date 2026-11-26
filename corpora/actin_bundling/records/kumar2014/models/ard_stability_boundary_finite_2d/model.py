import math


def run(inputs, settings):
    beta = inputs["beta"]
    K = inputs["kappa_tau"]
    rho = inputs["rho"]
    nu = inputs["nu"]
    L = inputs["L_over_l"]
    n_max = int(settings.get("n_max", 20))
    a = 1.0 + beta
    b = 1.0 - beta
    c_lin = K * (rho * a + b)
    dk = 2.0 * math.pi / L

    P_c_osc = float("inf")
    P_c_stat = float("inf")
    q_osc_best = None
    for nx in range(0, n_max + 1):
        for ny in range(0, n_max + 1):
            if nx == 0 and ny == 0:
                continue
            q = (dk**2) * (nx * nx + ny * ny)
            p_osc = (1.0 + nu * q) * (1.0 + K * (1.0 + rho) / (2.0 * q))
            if p_osc < P_c_osc:
                # verify det > 0 at this oscillatory threshold
                Pi = p_osc / (1.0 + nu * q)
                det = (q + K) * (q + rho * K) - q * Pi * (a * (q + rho * K) + b * (q + K))
                if det > 0.0:
                    P_c_osc = p_osc
                    q_osc_best = q
            denom = q * (2.0 * q + c_lin)
            if denom > 0.0:
                p_stat = (1.0 + nu * q) * (q + K) * (q + rho * K) / denom
                if p_stat < P_c_stat:
                    P_c_stat = p_stat

    if q_osc_best is not None and P_c_osc < P_c_stat:
        P_c = P_c_osc
        is_osc = 1.0
    else:
        P_c = P_c_stat
        is_osc = 0.0

    return {
        "P_c": P_c,
        "P_c_osc": P_c_osc,
        "P_c_stat": P_c_stat,
        "is_oscillatory": is_osc,
    }

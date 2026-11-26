import math
from scipy.optimize import minimize_scalar


def run(inputs, settings):
    beta = inputs["beta"]
    K = inputs["kappa_tau"]
    rho = inputs["rho"]
    nu = inputs["nu"]
    a = 1.0 + beta
    b = 1.0 - beta

    # Exact closed-form minimum of tr(tau*L_ARD) = 0 over q = k^2*ell^2 > 0
    q_osc = math.sqrt(K * (1.0 + rho) / (2.0 * nu))
    P_c_osc = (1.0 + math.sqrt(0.5 * nu * K * (1.0 + rho))) ** 2
    k_c_osc_l = math.sqrt(q_osc)

    # Exact minimum of det(tau*L_ARD) = 0 over q > q_min where denominator > 0
    # denom(q) = q * (2*q + K*(rho*a + b))
    c_lin = K * (rho * a + b)
    q_pole = max(0.0, -0.5 * c_lin)

    def p_stat_of_q(q):
        return (1.0 + nu * q) * (q + K) * (q + rho * K) / (q * (2.0 * q + c_lin))

    # Bracket minimum in log(q - q_pole)
    res = minimize_scalar(
        lambda u: p_stat_of_q(q_pole + math.exp(u)),
        bounds=(-12.0, 12.0),
        method="bounded",
        options={"xatol": 1e-12},
    )
    P_c_stat = float(res.fun)

    # Check whether oscillatory mode has det > 0 at (P_c_osc, q_osc)
    Pi_osc = P_c_osc / (1.0 + nu * q_osc)
    det_at_osc = (q_osc + K) * (q_osc + rho * K) - q_osc * Pi_osc * (
        a * (q_osc + rho * K) + b * (q_osc + K)
    )
    if det_at_osc > 0.0 and P_c_osc < P_c_stat:
        P_c = P_c_osc
        is_osc = 1.0
    else:
        P_c = P_c_stat
        is_osc = 0.0

    return {
        "P_c": P_c,
        "P_c_osc": P_c_osc,
        "P_c_stat": P_c_stat,
        "k_c_osc_l": k_c_osc_l,
        "is_oscillatory": is_osc,
    }

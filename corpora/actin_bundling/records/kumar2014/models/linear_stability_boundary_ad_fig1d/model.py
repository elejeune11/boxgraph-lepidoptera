import math


def run(inputs, settings):
    beta = inputs["beta"]
    alpha = inputs["alpha"]
    k1_l_sq = (inputs["omega"] * math.pi / inputs["L_over_l"]) ** 2
    fac = 1.0 + inputs["nu"] * k1_l_sq
    P_c_osc = 0.5 * (1.0 + alpha) * fac
    denom_stat = alpha * (1.0 + beta) + (1.0 - beta)
    if abs(1.0 - alpha) > 1e-12:
        beta_codim2 = (1.0 + alpha**2) / (1.0 - alpha**2)
    else:
        beta_codim2 = 1e9

    # Oscillatory precedes stationary iff det(L) > 0 at P = P_c_osc
    # i.e. alpha - ((1+alpha)/2) * (alpha*(1+beta) + (1-beta)) > 0
    det_at_osc = alpha - 0.5 * (1.0 + alpha) * denom_stat
    if det_at_osc > 0.0:
        P_c = P_c_osc
        is_osc = 1.0
    else:
        P_c = alpha * fac / denom_stat
        is_osc = 0.0
    return {
        "P_c": P_c,
        "P_c_osc": P_c_osc,
        "is_oscillatory": is_osc,
        "beta_codim2": beta_codim2,
    }

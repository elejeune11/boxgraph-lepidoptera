import math


def run(inputs, settings):
    Pe = inputs["Pe"]
    c0 = inputs["c0"]
    beta = inputs["beta"]
    L_over_l = inputs["L_over_l"]
    alpha = Pe * c0 / ((1.0 + c0) ** 2)
    dq = beta * math.pi / L_over_l

    if alpha <= 1.0:
        q1 = dq
        lam1 = -(q1**2) * (1.0 - alpha / (1.0 + q1**2))
        return {
            "n_max": 0.0,
            "k_max_discrete_l": 0.0,
            "lam_max_dimless": lam1,
            "n_unstable_modes": 0.0,
        }

    q_cut = math.sqrt(alpha - 1.0)
    n_upper = max(1, int(math.floor(q_cut / dq)) + 2)
    best_n = 1
    best_lam = -float("inf")
    n_unstable = 0
    for n in range(1, n_upper + 1):
        q = n * dq
        lam = -(q**2) * (1.0 - alpha / (1.0 + q**2))
        if lam > 0.0:
            n_unstable += 1
        if lam > best_lam:
            best_lam = lam
            best_n = n

    return {
        "n_max": float(best_n) if n_unstable > 0 else 0.0,
        "k_max_discrete_l": float(best_n * dq) if n_unstable > 0 else 0.0,
        "lam_max_dimless": float(best_lam),
        "n_unstable_modes": float(n_unstable),
    }

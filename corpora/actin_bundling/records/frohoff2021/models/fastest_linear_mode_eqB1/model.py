import math


def run(inputs, settings):
    a = inputs["a"]
    a_Delta = inputs["a_Delta"]
    kappa = inputs["kappa"]
    p1 = inputs["phi1_bar"]
    p2 = inputs["phi2_bar"]
    rho = inputs["rho"]
    alpha = inputs["alpha"]
    ell = inputs["ell"]
    n_max = int(settings.get("n_max", 25))

    Delta = alpha * alpha - rho * rho
    best_n = 1
    best_lam = -float("inf")
    bifs = []

    for n in range(1, n_max + 1):
        q = 2.0 * math.pi * n / ell
        q2 = q * q
        mean_term = -(q2 * (1.0 + kappa) + 2.0 * a + a_Delta + 3.0 * (p1 * p1 + p2 * p2))
        diff_term = q2 * (1.0 - kappa) + 3.0 * (p1 * p1 - p2 * p2) - a_Delta
        disc = diff_term * diff_term - 4.0 * Delta
        re_lam = 0.5 * q2 * (mean_term + math.sqrt(disc)) if disc >= 0.0 else 0.5 * q2 * mean_term
        if re_lam > best_lam:
            best_lam = re_lam
            best_n = n
        if disc >= 0.0:
            a_p = -0.5 * (q2 * (1.0 + kappa) + a_Delta + 3.0 * (p1 * p1 + p2 * p2)) + 0.5 * math.sqrt(disc)
            bifs.append((a_p, n))

    bifs.sort(key=lambda t: t[0], reverse=True)
    a_first, n_first = bifs[0] if len(bifs) >= 1 else (None, None)
    a_second, n_second = bifs[1] if len(bifs) >= 2 else (None, None)

    return {
        "n_fastest": float(best_n),
        "lam_fastest": best_lam,
        "n_first_bif": float(n_first) if n_first is not None else None,
        "a_first_bif": a_first,
        "n_second_bif": float(n_second) if n_second is not None else None,
        "a_second_bif": a_second,
    }

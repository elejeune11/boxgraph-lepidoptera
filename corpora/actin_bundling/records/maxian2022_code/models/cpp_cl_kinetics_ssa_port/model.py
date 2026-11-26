"""Faithful Python port of the 4-reaction event-driven Gillespie SSA in
Python/cppmodules/EndedCrossLinkers.cpp (lines 96-271, 524-539) coupled with
Python/DoubleEndedCrossLinkedNetwork.py (lines 96-121) on a fixed two-site pair
separated by distance r_sep."""

import math
import numpy as np


def _exact_gibbs_occupancy(lam1, lam2, max_per_site):
    """Exact stationary mean <m> and <n_0> of the 2-site Markov chain with
    P(n0, n1, m) proportional to lam1^(n0 + n1 + m) * lam2^m / (n0! n1! m!)
    restricted to n0 + m <= M and n1 + m <= M."""
    Z = 0.0
    m_sum = 0.0
    n0_sum = 0.0
    for m in range(max_per_site + 1):
        rem = max_per_site - m
        for n0 in range(rem + 1):
            for n1 in range(rem + 1):
                w = (
                    (lam1 ** (n0 + n1 + m))
                    * (lam2 ** m)
                    / (math.factorial(n0) * math.factorial(n1) * math.factorial(m))
                )
                Z += w
                m_sum += m * w
                n0_sum += n0 * w
    return n0_sum / Z, m_sum / Z


def run(inputs, settings):
    k_on = float(inputs["k_on"])
    k_off = float(inputs["k_off"])
    k_on_s = float(inputs["k_on_s"])
    k_off_s = float(inputs["k_off_s"])
    ds = float(inputs["ds"])
    c_w = float(inputs["c_w"])
    K_c = float(inputs["K_c"])
    ell = float(inputs["ell"])
    kBT = float(inputs["kBT"])
    r_sep = float(inputs["r_sep"])

    t_sim = float(settings.get("t_sim", 2000.0))
    burn_t = float(settings.get("burn_t", 50.0))
    seed = int(settings.get("seed", 0))

    # DoubleEndedCrossLinkedNetwork.py:96-99, 108-119
    kon_site = k_on * ds
    kons_site = k_on_s * ds
    max_per_site = int(math.ceil(ds / c_w)) if c_w > 0 else 1024

    if kBT > 0.0 and K_c > 0.0:
        deltaL = 2.0 * math.sqrt(kBT / K_c)
        energy = 0.5 * K_c * (r_sep - ell) ** 2
        boltz = math.exp(-energy / kBT)
    else:
        deltaL = min(math.sqrt(4.1e-3 / K_c), 0.5 * ell)
        boltz = 1.0

    lo = max(0.0, ell - deltaL) if deltaL <= ell else 0.0
    hi = ell + deltaL if deltaL <= ell else 2.0 * ell
    in_window = (r_sep > lo) and (r_sep < hi)

    # EndedCrossLinkers.cpp:524-530 + SpatialDatabase neighbor filter (DoubleEndedCrossLinkedNetwork.py:204)
    base_rate = (kons_site * boltz) if in_window else 0.0

    # Exact Gibbs-Boltzmann expectation WITHOUT capture-window truncation (Eq. 16-20)
    lam1 = kon_site / k_off
    lam2_exact = (kons_site * boltz) / k_off_s
    _, exact_m = _exact_gibbs_occupancy(lam1, lam2_exact, max_per_site)

    # Exact Gillespie SSA matching EndedCrossLinkers.cpp:125-270 for 2 sites (0 and 1)
    rng = np.random.default_rng(seed)
    free = [0, 0]
    m = 0  # doubly-bound links between site 0 and site 1
    t = 0.0
    acc_free = 0.0
    acc_m = 0.0
    obs_t = 0.0

    # Pre-bind constants
    r_free_bind = 2.0 * kon_site  # _kon * _TotSites (line 147)

    while t < t_sim:
        r_f0 = k_off * free[0]
        r_f1 = k_off * free[1]
        r_b01 = base_rate * free[0]
        r_b10 = base_rate * free[1]
        r_unb = 2.0 * k_off_s * m  # line 163
        r_tot = r_free_bind + r_f0 + r_f1 + r_b01 + r_b10 + r_unb
        if r_tot <= 0.0:
            break

        dt_ev = -math.log(1.0 - float(rng.random())) / r_tot
        t_next = min(t + dt_ev, t_sim)
        if t_next > burn_t:
            dt_obs = t_next - max(t, burn_t)
            acc_free += 0.5 * (free[0] + free[1]) * dt_obs
            acc_m += m * dt_obs
            obs_t += dt_obs
        t = t_next
        if t >= t_sim:
            break

        u = float(rng.random()) * r_tot
        if u < r_free_bind:
            # Reaction 1 (lines 172-180): random site, reject if TotalNumberBound >= max_per_site
            site = 0 if rng.random() < 0.5 else 1
            if free[site] + m < max_per_site:
                free[site] += 1
        elif u < r_free_bind + r_f0:
            # Reaction 2 at site 0 (lines 181-185)
            free[0] -= 1
        elif u < r_free_bind + r_f0 + r_f1:
            # Reaction 2 at site 1
            free[1] -= 1
        elif u < r_free_bind + r_f0 + r_f1 + r_b01:
            # Reaction 3: 0 -> 1 (lines 203-231)
            if free[1] + m < max_per_site:
                free[0] -= 1
                m += 1
        elif u < r_free_bind + r_f0 + r_f1 + r_b01 + r_b10:
            # Reaction 3: 1 -> 0
            if free[0] + m < max_per_site:
                free[1] -= 1
                m += 1
        else:
            # Reaction 4 (lines 186-203): unbind one random end, leaving free end on the other
            m -= 1
            if rng.random() < 0.5:
                free[0] += 1
            else:
                free[1] += 1

    mean_free = acc_free / obs_t if obs_t > 0 else 0.0
    mean_m = acc_m / obs_t if obs_t > 0 else 0.0
    rel_err = abs(mean_m - exact_m) / exact_m if exact_m > 0 else 0.0

    return {
        "site_binding_rate_at_r": float(base_rate),
        "mean_free_ends": float(mean_free),
        "mean_double_links": float(mean_m),
        "exact_gibbs_double_links": float(exact_m),
        "db_rel_error": float(rel_err),
    }

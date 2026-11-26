import math
import numpy as np
from scipy.optimize import curve_fit


def _simulate_runs(r, eps_kBT, n_sites, n_sims, seed):
    """Exact holding-time kinetic Monte Carlo on a 1D lattice of n_sites.

    Equivalent in distribution to site-by-site KMC (SI pp. 27-28):
    - First two sites are initialized identically as A (fascin) with probability
      r / (r + 1) and B (alpha-actinin) with probability 1 / (r + 1).
    - Each subsequent site switches from A to B with probability
      p_AB = 1 / (1 + r * exp(eps_kBT)) and from B to A with probability
      p_BA = r / (r + exp(eps_kBT)).
    """
    rng = np.random.default_rng(int(seed))
    N = int(round(n_sites))
    M = int(round(n_sims))
    p_A_init = r / (r + 1.0)
    exp_eps = math.exp(eps_kBT)
    p_AB = 1.0 / (1.0 + r * exp_eps)
    p_BA = r / (r + exp_eps)

    state = (rng.random(M) >= p_A_init).astype(np.int8)  # 0 = A (fascin), 1 = B (alpha-actinin)
    rem = np.full(M, N, dtype=np.int64)
    runs_A_list = []
    runs_B_list = []

    active = np.arange(M)
    first = True
    while active.size > 0:
        st = state[active]
        rm = rem[active]
        p_sw = np.where(st == 0, p_AB, p_BA)
        g = rng.geometric(p_sw)
        length = np.minimum(rm, g + 1 if first else g) if N >= 2 else rm
        first = False
        is_A = st == 0
        if np.any(is_A):
            runs_A_list.append(length[is_A])
        if np.any(~is_A):
            runs_B_list.append(length[~is_A])
        rem[active] -= length
        state[active] = 1 - st
        active = active[rem[active] > 0]

    runs_A = np.concatenate(runs_A_list) if runs_A_list else np.array([], dtype=np.int64)
    runs_B = np.concatenate(runs_B_list) if runs_B_list else np.array([], dtype=np.int64)
    return runs_A, runs_B


def _fit_exp_decay(runs, n_sites, bin_width_sites, a_um):
    bw = max(1, int(round(bin_width_sites)))
    N = int(round(n_sites))
    edges = np.arange(1, N + bw + 1, bw)
    if edges[-1] <= N:
        edges = np.append(edges, edges[-1] + bw)
    centers_um = (0.5 * (edges[:-1] + edges[1:]) - 0.5) * a_um
    hist, _ = np.histogram(runs, bins=edges, density=True)
    L0 = max(float(runs.mean()) * a_um, a_um)
    A0 = max(float(hist[0]), 1e-9)
    popt, _ = curve_fit(
        lambda x, A, L: A * np.exp(-x / L),
        centers_um,
        hist,
        p0=[A0, L0],
        bounds=([0.0, 1e-6], [np.inf, np.inf]),
        maxfev=20000,
    )
    return float(popt[1])


def run(inputs, settings):
    r = float(inputs["r"])
    eps_kBT = float(inputs["eps_kBT"])
    a_um = float(inputs["a"]) * 1e-3
    n_sites = int(round(settings.get("n_sites", 400)))
    n_sims = int(round(settings.get("n_sims", 5000)))
    bw = int(round(settings.get("bin_width_sites", 10)))
    seed = int(settings.get("seed", 0))

    if r <= 0.0:
        return {"l_F": 0.0, "l_alpha": float(n_sites) * a_um}

    runs_A, runs_B = _simulate_runs(r, eps_kBT, n_sites, n_sims, seed)
    return {
        "l_F": _fit_exp_decay(runs_A, n_sites, bw, a_um),
        "l_alpha": _fit_exp_decay(runs_B, n_sites, bw, a_um),
    }

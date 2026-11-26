import math
import numpy as np


def _simulate_runs(r, eps_kBT, n_sites, n_sims, seed):
    """Exact holding-time kinetic Monte Carlo on a 1D lattice of n_sites (SI pp. 27-28)."""
    rng = np.random.default_rng(int(seed))
    N = int(round(n_sites))
    M = int(round(n_sims))
    p_A_init = r / (r + 1.0)
    exp_eps = math.exp(eps_kBT)
    p_AB = 1.0 / (1.0 + r * exp_eps)
    p_BA = r / (r + exp_eps)

    state = (rng.random(M) >= p_A_init).astype(np.int8)
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


def run(inputs, settings):
    r = float(inputs["r"])
    eps_kBT = float(inputs["eps_kBT"])
    a_um = float(inputs["a"]) * 1e-3
    n_sites = int(round(settings.get("n_sites", 400)))
    n_sims = int(round(settings.get("n_sims", 5000)))
    seed = int(settings.get("seed", 0))

    if r <= 0.0:
        return {"l_F": 0.0, "l_alpha": float(n_sites) * a_um}

    runs_A, runs_B = _simulate_runs(r, eps_kBT, n_sites, n_sims, seed)
    return {
        "l_F": float(runs_A.mean()) * a_um,
        "l_alpha": float(runs_B.mean()) * a_um,
    }

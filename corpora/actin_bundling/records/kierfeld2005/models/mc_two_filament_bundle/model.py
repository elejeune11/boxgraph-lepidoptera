"""Metropolis Monte Carlo simulation of two identical parallel semiflexible filaments
under the effective Hamiltonian H = H_{b,1} + H_{b,2} + H_bar_2 (Kierfeld et al., PRL 2005).

Uses the exact decoupling of the center-of-mass mode (p. 2 col 1: mean CM bending energy
(L - 1) T in 2D) and simulates the relative displacement field r_k = z_{1,k} - z_{2,k} with
reduced rigidity kappa_{12}^r = kappa / 2 = L_p T / 4 (from L_p = 2 kappa / T, p. 1 col 2),
periodic boundary conditions along x, 4-sublattice single-bead updates, and collective
bump / global translation moves obeying detailed balance.
"""
import numpy as np


def simulate_pair(W_abs, L_p, L, l_a, l_r, n_equil, n_meas, seed):
    rng = np.random.default_rng(int(seed))
    L = int(round(L))
    n_equil = int(round(n_equil))
    n_meas = int(round(n_meas))
    kappa_r = 0.25 * float(L_p)

    # Start bound inside the attractive square well
    r = np.zeros((L, 2), dtype=float)
    r[:, 0] = l_r + 0.5 * l_a

    sub = [np.arange(s, L, 4) for s in range(4)]
    m2 = [(s - 2) % L for s in sub]
    m1 = [(s - 1) % L for s in sub]
    p1 = [(s + 1) % L for s in sub]
    p2 = [(s + 2) % L for s in sub]
    steps = (0.0015, 0.008, 0.05, 0.30)
    lr2 = l_r * l_r
    lra2 = (l_r + l_a) ** 2

    e_sum = 0.0
    d1_sum = 0.0
    d2_sum = 0.0
    n_rec = 0

    for sw in range(n_equil + n_meas):
        step = steps[sw & 3]
        for si in range(4):
            idx = sub[si]
            old = r[idx]
            prop = old + rng.uniform(-step, step, size=old.shape)
            ro2 = old[:, 0] * old[:, 0] + old[:, 1] * old[:, 1]
            rn2 = prop[:, 0] * prop[:, 0] + prop[:, 1] * prop[:, 1]
            ov = rn2 < lr2
            in_old = (ro2 >= lr2) & (ro2 < lra2)
            in_new = (rn2 >= lr2) & (rn2 < lra2)
            dV = -W_abs * (in_new.astype(float) - in_old.astype(float))
            zm2, zm1, zp1, zp2 = r[m2[si]], r[m1[si]], r[p1[si]], r[p2[si]]
            c1o = old - 2.0 * zm1 + zm2
            c0o = zp1 - 2.0 * old + zm1
            c2o = zp2 - 2.0 * zp1 + old
            c1n = prop - 2.0 * zm1 + zm2
            c0n = zp1 - 2.0 * prop + zm1
            c2n = zp2 - 2.0 * zp1 + prop
            dHb = 0.5 * kappa_r * np.sum(
                c1n * c1n + c0n * c0n + c2n * c2n - c1o * c1o - c0o * c0o - c2o * c2o,
                axis=1,
            )
            dE = dHb + dV
            acc = (~ov) & ((dE <= 0.0) | (rng.random(len(idx)) < np.exp(-np.clip(dE, 0.0, 500.0))))
            r[idx] = np.where(acc[:, None], prop, old)

        if (sw & 3) == 0:
            k0 = int(rng.integers(0, L))
            w = int(rng.integers(4, L))
            if w > L - 6:
                amp = rng.normal(0.0, 0.25, size=(1, 2))
                r_new = r + amp
                rn2 = r_new[:, 0] ** 2 + r_new[:, 1] ** 2
                if not np.any(rn2 < lr2):
                    ro2 = r[:, 0] ** 2 + r[:, 1] ** 2
                    dE = -W_abs * (
                        np.sum((rn2 >= lr2) & (rn2 < lra2))
                        - np.sum((ro2 >= lr2) & (ro2 < lra2))
                    )
                    if dE <= 0.0 or rng.random() < np.exp(-min(dE, 500.0)):
                        r = r_new
            else:
                amp = rng.normal(0.0, 0.015 * (w / 10.0) ** 1.5, size=(1, 2))
                ks = (k0 + np.arange(w)) % L
                prof = np.sin(np.pi * (np.arange(w) + 0.5) / w)[:, None]
                r_new = r.copy()
                r_new[ks] += amp * prof
                rn2 = r_new[ks, 0] ** 2 + r_new[ks, 1] ** 2
                if not np.any(rn2 < lr2):
                    ro2 = r[ks, 0] ** 2 + r[ks, 1] ** 2
                    dV = -W_abs * (
                        np.sum((rn2 >= lr2) & (rn2 < lra2))
                        - np.sum((ro2 >= lr2) & (ro2 < lra2))
                    )
                    ck = (k0 - 1 + np.arange(min(w + 2, L))) % L
                    cm = (ck - 1) % L
                    cp = (ck + 1) % L
                    co = r[cp] - 2.0 * r[ck] + r[cm]
                    cn = r_new[cp] - 2.0 * r_new[ck] + r_new[cm]
                    dE = 0.5 * kappa_r * (np.sum(cn * cn) - np.sum(co * co)) + dV
                    if dE <= 0.0 or rng.random() < np.exp(-min(dE, 500.0)):
                        r = r_new

        if sw >= n_equil and (sw % 20) == 0:
            c = np.roll(r, -1, axis=0) - 2.0 * r + np.roll(r, 1, axis=0)
            hb_rel = 0.5 * kappa_r * np.sum(c * c)
            rn = np.sqrt(r[:, 0] ** 2 + r[:, 1] ** 2)
            hv = -W_abs * np.sum((rn >= l_r) & (rn < l_r + l_a))
            # Center-of-mass mode contributes exactly (L - 1) T by 2D quadratic equipartition
            e_sum += ((L - 1.0) + hb_rel + hv) / (2.0 * L)
            dr = rn - l_r
            d1_sum += float(np.mean(dr))
            d2_sum += float(np.mean(dr * dr))
            n_rec += 1

    return e_sum / n_rec, d1_sum / n_rec, float(np.sqrt(d2_sum / n_rec))


def run(inputs, settings):
    e_mean, dz_mean, _ = simulate_pair(
        W_abs=inputs["W_bar"],
        L_p=inputs["L_p"],
        L=inputs["L"],
        l_a=inputs["l_a"],
        l_r=inputs["l_r"],
        n_equil=settings.get("n_equil", 3000),
        n_meas=settings.get("n_meas", 4000),
        seed=settings.get("seed", 1),
    )
    return {
        "energy_per_filament_length": float(e_mean),
        "mean_separation": float(dz_mean),
    }

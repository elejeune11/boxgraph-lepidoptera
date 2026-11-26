import math
import random


def _extract_domains(state, lb_nm, n_d, ld_nm):
    short_doms = []
    long_doms = []
    cur_t = 0
    first_pos = 0
    last_pos = 0
    count = 0
    for idx, v in enumerate(state):
        if v == 0:
            continue
        if cur_t == 0:
            cur_t = v
            first_pos = idx
            last_pos = idx
            count = 1
        elif v == cur_t and (idx - last_pos) * lb_nm < ld_nm:
            last_pos = idx
            count += 1
        else:
            if count >= n_d:
                length_um = (last_pos - first_pos) * lb_nm * 1e-3
                if cur_t == 1:
                    short_doms.append(length_um)
                else:
                    long_doms.append(length_um)
            cur_t = v
            first_pos = idx
            last_pos = idx
            count = 1
    if cur_t != 0 and count >= n_d:
        length_um = (last_pos - first_pos) * lb_nm * 1e-3
        if cur_t == 1:
            short_doms.append(length_um)
        else:
            long_doms.append(length_um)
    return short_doms, long_doms


def run(inputs, settings):
    N = int(round((inputs["L_0"] * 1000.0) / inputs["l_b"]))
    lb_nm = float(inputs["l_b"])
    Lp_nm = float(inputs["L_p"]) * 1000.0
    dl_nm = float(inputs["delta_l_xl"])
    mu_S = float(inputs["mu_S_kBT"])
    mu_L = float(inputs["mu_L_kBT"])
    n_d = int(round(inputs["n_d"]))
    ld_nm = float(inputs["l_d"]) * 1000.0

    n_sweeps = int(settings.get("n_sweeps", 200))
    n_sims = int(settings.get("n_sims", 20))
    seed = int(settings.get("seed", 0))

    # Precompute Eq. (1) gap energy U(d * lb) / (k_B T).
    # Gaps with contour length <= delta_l_xl are geometrically impossible (infinite energy).
    U_tab = [0.0] * (N + 1)
    for d in range(1, N + 1):
        lg = d * lb_nm
        if lg <= dl_nm:
            U_tab[d] = 1e9
        else:
            th = math.asin(dl_nm / lg)
            U_tab[d] = (Lp_nm / lg) * (th * th)

    mu_map = (0.0, mu_S, mu_L)
    rng = random.Random(seed)
    rand = rng.random
    exp = math.exp

    all_short = []
    all_long = []
    n_steps = n_sweeps * N

    for _ in range(n_sims):
        state = [0] * N
        for _ in range(n_steps):
            i = int(rand() * N)
            old_s = state[i]
            new_s = (old_s + (1 if rand() < 0.5 else 2)) % 3

            lp = i - 1
            while lp >= 0 and state[lp] == 0:
                lp -= 1
            rp = i + 1
            while rp < N and state[rp] == 0:
                rp += 1

            lt = state[lp] if lp >= 0 else 0
            rt = state[rp] if rp < N else 0

            u_old = 0.0
            if old_s == 0:
                if lt and rt and lt != rt:
                    u_old = U_tab[rp - lp]
            else:
                if lt and lt != old_s:
                    u_old += U_tab[i - lp]
                if rt and rt != old_s:
                    u_old += U_tab[rp - i]

            u_new = 0.0
            if new_s == 0:
                if lt and rt and lt != rt:
                    u_new = U_tab[rp - lp]
            else:
                if lt and lt != new_s:
                    u_new += U_tab[i - lp]
                if rt and rt != new_s:
                    u_new += U_tab[rp - i]

            dU = -(mu_map[new_s] - mu_map[old_s]) + (u_new - u_old)
            if dU <= 0.0 or (dU < 700.0 and rand() < exp(-dU)):
                state[i] = new_s

        s_doms, l_doms = _extract_domains(state, lb_nm, n_d, ld_nm)
        all_short.extend(s_doms)
        all_long.extend(l_doms)

    all_doms = all_short + all_long
    mean_all = sum(all_doms) / len(all_doms) if all_doms else 0.0
    mean_s = sum(all_short) / len(all_short) if all_short else 0.0
    mean_l = sum(all_long) / len(all_long) if all_long else 0.0
    return {
        "domain_length": mean_all,
        "domain_length_short": mean_s,
        "domain_length_long": mean_l,
    }

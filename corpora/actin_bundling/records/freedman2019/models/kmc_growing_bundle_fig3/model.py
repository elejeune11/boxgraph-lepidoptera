import math
import random


def run(inputs, settings):
    mu_grow = float(inputs["mu_grow"])
    mu_A = float(inputs["mu_on_A"])
    mu_B = float(inputs["mu_on_B"])
    lb_nm = float(inputs["l_b"])
    Lp_nm = float(inputs["L_p"]) * 1000.0
    dl_nm = float(inputs["delta_l_xl"])
    N_max = int(round((inputs["L_0"] * 1000.0) / lb_nm))
    n_d = int(round(inputs["n_d"]))
    ld_nm = float(inputs["l_d"]) * 1000.0

    n_sims = int(settings.get("n_sims", 1000))
    seed = int(settings.get("seed", 0))

    r_grow = mu_grow / lb_nm
    expU = [0.0] * (N_max + 2)
    for d in range(1, N_max + 2):
        lg = d * lb_nm
        if lg > dl_nm:
            th = math.asin(dl_nm / lg)
            u = (Lp_nm / lg) * (th * th)
            expU[d] = math.exp(-u) if u < 700.0 else 0.0

    # Prefix sums of expU[1..m] for O(1) total dissimilar binding rate
    expU_cum = [0.0] * (N_max + 2)
    acc = 0.0
    for d in range(1, N_max + 2):
        acc += expU[d]
        expU_cum[d] = acc

    rng = random.Random(seed)
    rand = rng.random
    p_init_A = mu_A / (mu_A + mu_B)

    n_same = 0
    n_diff = 0
    dom_A = []
    dom_B = []

    for _ in range(n_sims):
        cur_t = 1 if rand() < p_init_A else 2
        last_b = 0
        cur_len = 1
        bound = [(0, cur_t)]

        while cur_len < N_max:
            m = cur_len - 1 - last_b
            mu_s = mu_A if cur_t == 1 else mu_B
            mu_d = mu_B if cur_t == 1 else mu_A
            r_same = m * mu_s
            r_diff = mu_d * expU_cum[m]
            R_tot = r_grow + r_same + r_diff
            u = rand() * R_tot

            if u < r_grow:
                cur_len += 1
            elif u < r_grow + r_same:
                d = int((u - r_grow) / mu_s) + 1
                if d > m:
                    d = m
                last_b += d
                bound.append((last_b, cur_t))
                n_same += 1
            else:
                target = (u - r_grow - r_same) / mu_d
                # Linear scan over small m (typically m < 20)
                d = 1
                while d < m and expU_cum[d] < target:
                    d += 1
                last_b += d
                cur_t = 2 if cur_t == 1 else 1
                bound.append((last_b, cur_t))
                n_diff += 1

        # SI Section B domain extraction: consecutive run of >= n_d cross-linkers separated by < l_d
        ct = bound[0][1]
        first_pos = bound[0][0]
        last_pos = first_pos
        count = 1
        for pos, tp in bound[1:]:
            if tp == ct and (pos - last_pos) * lb_nm < ld_nm:
                last_pos = pos
                count += 1
            else:
                if count >= n_d:
                    L_dom = (last_pos - first_pos) * lb_nm * 1e-3
                    if ct == 1:
                        dom_A.append(L_dom)
                    else:
                        dom_B.append(L_dom)
                ct = tp
                first_pos = pos
                last_pos = pos
                count = 1
        if count >= n_d:
            L_dom = (last_pos - first_pos) * lb_nm * 1e-3
            if ct == 1:
                dom_A.append(L_dom)
            else:
                dom_B.append(L_dom)

    ratio = n_diff / n_same if n_same > 0 else float("nan")
    eps_pred = -math.log(ratio) if ratio > 0.0 else float("inf")
    mean_A = sum(dom_A) / len(dom_A) if dom_A else 0.0
    mean_B = sum(dom_B) / len(dom_B) if dom_B else 0.0
    return {
        "switching_ratio": ratio,
        "epsilon_predicted": eps_pred,
        "domain_length_A": mean_A,
        "domain_length_B": mean_B,
    }

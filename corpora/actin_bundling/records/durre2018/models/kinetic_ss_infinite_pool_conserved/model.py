"""Exact 1D infinite-pool steady state of Eqs. (6)-(16) with mass-conserving boundary conditions (+k_2*V* in Eq. 14 and +k_4*A* in Eq. 15)."""

from scipy.optimize import brentq

NA_PER_UM3_PER_UM = 602.214076


def run(inputs, settings):
    G0 = float(inputs["G_0"])
    Arp0 = float(inputs["Arp_0"])
    CP0 = float(inputs["CP_0"]) * 1e-3  # nM -> uM
    V0 = float(inputs["V_0"])
    h = float(inputs["h"])
    kon = float(inputs["k_on"])
    koff = float(inputs["k_off"])
    k1 = float(inputs["k_1"])
    k2 = float(inputs["k_2"])
    k3 = float(inputs["k_3"]) / float(inputs["k_3_reduction"])
    k4 = float(inputs["k_4"])
    kCap = float(inputs["k_Cap"])
    kArp = float(inputs["k_Arp"])
    kre = float(inputs["k_re"])
    kDeg = float(inputs["k_Deg"])
    DG = float(inputs["D_G"])
    DArp = float(inputs["D_Arp"])
    DCP = float(inputs["D_CP"])
    xtol = float(settings.get("xtol", 1e-12))

    # Trivial B = 0 equilibrium (no net flux across membrane when reverse terms are included)
    K1 = k1 * G0 / k2
    K3 = k3 * Arp0 / k4
    V_eq0 = V0 / (1.0 + K1 + K1 * K3)
    As_eq0 = K1 * K3 * V_eq0

    if kArp * As_eq0 <= kCap * CP0 + kDeg:
        return {
            "B_ss": 0.0,
            "A_star_ss": float(As_eq0),
            "branching_rate_ss": 0.0,
            "G_mem_ss": float(G0),
            "Arp_mem_ss": float(Arp0),
            "CP_mem_ss": float(CP0 * 1e3),
        }

    def eval_at_B(B):
        CPm = CP0 / (1.0 + kCap * h * B / (DCP * NA_PER_UM3_PER_UM))
        As = (kCap * CPm + kDeg) / kArp
        br = kArp * As * B
        Vp = br / kre
        Arpm = Arp0 - br * h / (DArp * NA_PER_UM3_PER_UM)
        if Arpm <= 0.0:
            return None
        Vs = (k4 * As + br) / (k3 * Arpm)
        V = V0 - Vs - As - Vp
        if V <= 0.0:
            return None
        Gm = (k2 * Vs + br) / (k1 * V)
        flux_G = DG * (G0 - Gm) * NA_PER_UM3_PER_UM / h
        sink_G = br + (kon * Gm - koff) * B
        return flux_G - sink_G, (As, br, Gm, Arpm, CPm)

    b_lo = 1e-9
    b_hi = 1.0
    while True:
        out = eval_at_B(b_hi)
        if out is None or out[0] <= 0.0:
            break
        b_lo = b_hi
        b_hi *= 2.0
        if b_hi > 1e8:
            break

    if eval_at_B(b_hi) is None:
        lo_f, hi_f = b_lo, b_hi
        for _ in range(80):
            mid = 0.5 * (lo_f + hi_f)
            out = eval_at_B(mid)
            if out is None:
                hi_f = mid
            elif out[0] > 0.0:
                lo_f = mid
            else:
                b_lo, b_hi = lo_f, mid
                break
        else:
            b_lo, b_hi = lo_f, lo_f

    B_ss = brentq(lambda b: eval_at_B(b)[0], b_lo, b_hi, xtol=xtol) if b_hi > b_lo else b_lo
    _, (As, br, Gm, Arpm, CPm) = eval_at_B(B_ss)
    return {
        "B_ss": float(B_ss),
        "A_star_ss": float(As),
        "branching_rate_ss": float(br),
        "G_mem_ss": float(Gm),
        "Arp_mem_ss": float(Arpm),
        "CP_mem_ss": float(CPm * 1e3),
    }

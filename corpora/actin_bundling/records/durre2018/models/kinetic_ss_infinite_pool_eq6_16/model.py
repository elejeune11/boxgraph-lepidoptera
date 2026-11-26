"""Exact 1D infinite-pool steady state of printed Eqs. (6)-(16) (SI p. 12).

Implements printed Eqs. (14)-(15) verbatim (without reverse dissociation release
terms +k_2*V* and +k_4*A* in the cytosolic flux boundary conditions).
"""

from scipy.optimize import brentq

NA_PER_UM3_PER_UM = 602.214076


def _solve_steady_state(inputs, settings, conserved_bc):
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

    # 1) Solve trivial B = 0 state first to test linear growth instability
    if conserved_bc:
        Gm0, Arpm0, CPm0 = G0, Arp0, CP0
        K1 = k1 * Gm0 / k2
        K3 = k3 * Arpm0 / k4
        V_eq0 = V0 / (1.0 + K1 + K1 * K3)
        Vs_eq0 = K1 * V_eq0
        As_eq0 = K1 * K3 * V_eq0
    else:
        # Printed Eqs. (14)-(15) drain G and Arp even at B = 0:
        # DG*(G0 - Gm)/h * NA = k1*Gm*V = k2*Vs => Vs = DG*(G0 - Gm)*NA/(k2*h), V = k2*Vs/(k1*Gm)
        # DArp*(Arp0 - Arpm)/h * NA = k3*Arpm*Vs => Arpm = Arp0 / (1 + k3*Vs*h/(DArp*NA))
        def res_b0(Gm):
            Vs = DG * (G0 - Gm) * NA_PER_UM3_PER_UM / (k2 * h)
            V = k2 * Vs / (k1 * Gm)
            Arpm = Arp0 / (1.0 + k3 * Vs * h / (DArp * NA_PER_UM3_PER_UM))
            As = (k3 / k4) * Arpm * Vs
            return V + Vs + As - V0

        Gm0 = brentq(res_b0, 1e-12 * G0, G0 * (1.0 - 1e-15), xtol=xtol)
        Vs_eq0 = DG * (G0 - Gm0) * NA_PER_UM3_PER_UM / (k2 * h)
        V_eq0 = k2 * Vs_eq0 / (k1 * Gm0)
        Arpm0 = Arp0 / (1.0 + k3 * Vs_eq0 * h / (DArp * NA_PER_UM3_PER_UM))
        As_eq0 = (k3 / k4) * Arpm0 * Vs_eq0
        CPm0 = CP0

    # Check if B = 0 is stable (net growth rate <= 0)
    if kArp * As_eq0 <= kCap * CPm0 + kDeg:
        return {
            "B_ss": 0.0,
            "A_star_ss": float(As_eq0),
            "branching_rate_ss": 0.0,
            "G_mem_ss": float(Gm0),
            "Arp_mem_ss": float(Arpm0),
            "CP_mem_ss": float(CPm0 * 1e3),
        }

    # 2) For B > 0, express all species explicitly as functions of B and solve Eq. (14)
    def eval_at_B(B):
        CPm = CP0 / (1.0 + kCap * h * B / (DCP * NA_PER_UM3_PER_UM))
        As = (kCap * CPm + kDeg) / kArp
        br = kArp * As * B
        Vp = br / kre
        arp_sink = br if conserved_bc else (k4 * As + br)
        Arpm = Arp0 - arp_sink * h / (DArp * NA_PER_UM3_PER_UM)
        if Arpm <= 0.0:
            return None
        Vs = (k4 * As + br) / (k3 * Arpm)
        V = V0 - Vs - As - Vp
        if V <= 0.0:
            return None
        Gm = (k2 * Vs + br) / (k1 * V)
        flux_G = DG * (G0 - Gm) * NA_PER_UM3_PER_UM / h
        sink_G = (br if conserved_bc else (k2 * Vs + br)) + (kon * Gm - koff) * B
        return flux_G - sink_G, (As, br, Gm, Arpm, CPm)

    # Find upper bracket where V(B) -> 0 or Arpm(B) -> 0 or residual becomes negative
    b_lo = 1e-9
    r_lo = eval_at_B(b_lo)
    if r_lo is None or r_lo[0] <= 0.0:
        return {
            "B_ss": 0.0,
            "A_star_ss": float(As_eq0),
            "branching_rate_ss": 0.0,
            "G_mem_ss": float(Gm0),
            "Arp_mem_ss": float(Arpm0),
            "CP_mem_ss": float(CPm0 * 1e3),
        }

    b_hi = 1.0
    while True:
        out = eval_at_B(b_hi)
        if out is None or out[0] <= 0.0:
            break
        b_lo = b_hi
        b_hi *= 2.0
        if b_hi > 1e8:
            break

    # Tighten upper bracket if eval_at_B(b_hi) is infeasible (None)
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

    if b_hi > b_lo:
        B_ss = brentq(lambda b: eval_at_B(b)[0], b_lo, b_hi, xtol=xtol)
    else:
        B_ss = b_lo
    _, (As, br, Gm, Arpm, CPm) = eval_at_B(B_ss)
    return {
        "B_ss": float(B_ss),
        "A_star_ss": float(As),
        "branching_rate_ss": float(br),
        "G_mem_ss": float(Gm),
        "Arp_mem_ss": float(Arpm),
        "CP_mem_ss": float(CPm * 1e3),
    }


def run(inputs, settings):
    return _solve_steady_state(inputs, settings, conserved_bc=False)

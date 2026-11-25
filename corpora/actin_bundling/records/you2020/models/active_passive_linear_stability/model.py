import math


def run(inputs, settings):
    v0 = float(inputs["v_0"])
    DR = float(inputs["D_R"])
    a0AA = float(inputs["a0_AA"])
    a0AP = float(inputs["a0_AP"])
    a1 = float(inputs["a1"])
    kap = float(inputs["kappa_ap"])
    rA = float(inputs["rho_A0"])
    rP = float(inputs["rho_P0"])
    L = float(inputs["L"])

    v = v0 * (1.0 - a0AA * rA - a0AP * rP)
    DAA = a1 * rA + (v / (2.0 * DR)) * (v - a0AA * v0 * rA)
    DAP = a1 * rA - (v / (2.0 * DR)) * a0AP * v0 * rA
    DPA = a1 * rP
    DPP = a1 * rP

    # Eigenvalues of -[[DAA, DAP], [DPA, DPP]]:
    tr = -(DAA + DPP)
    diff = -(DAA - DPP)
    disc = diff * diff + 4.0 * DAP * DPA
    if disc >= 0.0:
        max_re_negD = 0.5 * (tr + math.sqrt(disc))
    else:
        max_re_negD = 0.5 * tr

    q1 = 2.0 * math.pi / L
    q1_2 = q1 * q1
    max_re_lam = q1_2 * (max_re_negD - kap * q1_2)
    L_min = 2.0 * math.pi * math.sqrt(kap / max_re_negD) if max_re_negD > 0.0 else None

    return {
        "D_AA": DAA,
        "D_AP": DAP,
        "D_PA": DPA,
        "D_PP": DPP,
        "max_re_lambda_q1": max_re_lam,
        "L_min_unstable": L_min,
    }

import math


def run(inputs, settings):
    aA = float(inputs["alpha_A"])
    aB = float(inputs["alpha_B"])
    k = float(inputs["kappa"])
    d = float(inputs["delta_over_kappa"]) * k
    dc = math.sqrt(k * k + aB * aB)

    arg_s = (k * k - d * d - aA * aB) / aB
    if arg_s >= 0.0:
        rho_A_s = math.sqrt(arg_s)
        rho_B_s = (k + d) * rho_A_s / aB
    else:
        rho_A_s = None
        rho_B_s = None

    if d < dc - 1e-14 or (-aA - aB) < 0.0:
        v = 0.0
        rho_A_t = None
        rho_B_t = None
        psi_t = 0.0
    else:
        v = math.sqrt(max(0.0, d * d - dc * dc))
        rho_A_t = math.sqrt(-aA - aB)
        rho_B_t = rho_A_t * math.sqrt((d + k) / (d - k))
        cos_arg = -min(1.0, math.sqrt(aB * aB / (d * d - k * k)))
        psi_t = math.pi - math.acos(cos_arg)

    return {
        "rho_A_s": rho_A_s,
        "rho_B_s": rho_B_s,
        "rho_A_t": rho_A_t,
        "rho_B_t": rho_B_t,
        "v": v,
        "psi_t": psi_t,
    }

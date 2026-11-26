import math

KB_PN_UM_PER_K = 1.380649e-5


def run(inputs, settings):
    iL = 1 if inputs["I_L"] >= 0.5 else 0
    iR = 1 if inputs["I_R"] >= 0.5 else 0
    u_stretch = 0.5 * inputs["k_xl"] * ((inputs["r"] - inputs["l_xl"]) ** 2)
    half_pi = 0.5 * math.pi
    u_bend = inputs["kappa_xl_B"] * (
        ((inputs["alpha_L"] - half_pi) ** 2) * iL
        + ((inputs["alpha_R"] - half_pi) ** 2) * iR
    )
    kBT = KB_PN_UM_PER_K * inputs["T"]
    any_bound = 1 if (iL or iR) else 0
    both_bound = 1 if (iL and iR) else 0
    u_bind = -kBT * (
        math.log(inputs["k_on_1"] / inputs["k_off_1"]) * any_bound
        + math.log(inputs["k_on_2"] / inputs["k_off_2"]) * both_bound
    )
    return {
        "U_stretch": u_stretch,
        "U_bend": u_bend,
        "U_bind": u_bind,
        "U_xl": u_stretch + u_bend + u_bind,
    }

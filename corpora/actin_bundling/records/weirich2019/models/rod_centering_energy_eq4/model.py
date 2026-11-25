import math


def _eq4_dimless(x_c, L_D, L_M):
    if x_c <= -0.5 * L_M:
        return -math.log(L_D ** 2 - (-0.5 * L_M + x_c) ** 2) + math.log(L_D ** 2 - (0.5 * L_M + x_c) ** 2)
    elif x_c <= 0.5 * L_M:
        return (
            math.log(L_D ** 4)
            - math.log(L_D ** 2 - ((L_M - 2.0 * x_c) ** 2) / 4.0)
            - math.log(L_D ** 2 - ((L_M + 2.0 * x_c) ** 2) / 4.0)
        )
    else:
        return (
            math.log(L_D ** 2 - ((L_M - 2.0 * x_c) ** 2) / 4.0)
            - math.log(L_D ** 2 - ((L_M + 2.0 * x_c) ** 2) / 4.0)
        )


def run(inputs, settings):
    F_0 = inputs["F_0"]  # pN*nm
    x_c = inputs["x_c"]  # um
    L_D = inputs["L_D"]  # um
    L_M = inputs["L_M"]  # um

    dimless = _eq4_dimless(x_c, L_D, L_M)
    dimless_0 = _eq4_dimless(0.0, L_D, L_M)
    shifted_dimless = dimless - dimless_0
    return {
        "Delta_F": F_0 * dimless,
        "Delta_F_over_F0": dimless,
        "Delta_F_shifted": F_0 * shifted_dimless,
        "Delta_F_shifted_over_F0": shifted_dimless,
    }

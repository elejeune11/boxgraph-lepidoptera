def run(inputs, settings):
    F_0_pN_um = inputs["F_0"] * 1e-3  # pN*nm -> pN*um
    x_abs = abs(inputs["x_c"])        # um (symmetric potential)
    L_D = inputs["L_D"]               # um
    L_M = inputs["L_M"]               # um
    Gamma = inputs["Gamma"]           # pN*s/um

    # Exact derivative of Eq. 4 for x_c >= 0 (continuous at x_c = L_M / 2)
    if x_abs <= 0.5 * L_M:
        dF_dimless_dx = (
            -2.0 * (0.5 * L_M - x_abs) / (L_D ** 2 - (0.5 * L_M - x_abs) ** 2)
            + 2.0 * (0.5 * L_M + x_abs) / (L_D ** 2 - (0.5 * L_M + x_abs) ** 2)
        )
    else:
        dF_dimless_dx = (
            -2.0 * (x_abs - 0.5 * L_M) / (L_D ** 2 - (x_abs - 0.5 * L_M) ** 2)
            + 2.0 * (x_abs + 0.5 * L_M) / (L_D ** 2 - (x_abs + 0.5 * L_M) ** 2)
        )

    f_pN = F_0_pN_um * dF_dimless_dx
    v_nm_s = (f_pN / Gamma) * 1000.0  # um/s -> nm/s
    return {"f": f_pN, "v": v_nm_s}

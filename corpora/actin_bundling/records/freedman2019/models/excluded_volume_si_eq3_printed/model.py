def run(inputs, settings):
    r_ij = inputs["r_ij"]
    r_max = inputs["r_max"]
    if r_ij < r_max:
        u = inputs["k_exv"] * ((1.0 - r_ij / r_max) ** 2)
    else:
        u = 0.0
    return {"U_ij": u}

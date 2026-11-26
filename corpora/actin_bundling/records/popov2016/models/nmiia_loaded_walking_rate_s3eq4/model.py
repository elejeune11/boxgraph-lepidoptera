def run(inputs, settings):
    f_per_head = inputs["F_ext"] / inputs["N_total"]
    k_w = inputs["k_ens_walk_0"] * (inputs["F_stall"] - f_per_head) / (
        inputs["F_stall"] + f_per_head / inputs["zeta"]
    )
    return {"k_ens_walk": k_w, "v_walk": k_w * inputs["d_step"]}

def symbolic():
    import sympy as sp
    k_ens_walk_0, F_stall, F_ext, N_total, zeta, d_step = sp.symbols(
        "k_ens_walk_0 F_stall F_ext N_total zeta d_step", positive=True
    )
    f_per_head = F_ext / N_total
    k_w = k_ens_walk_0 * (F_stall - f_per_head) / (F_stall + f_per_head / zeta)
    return {"k_ens_walk": k_w, "v_walk": k_w * d_step}

def run(inputs, settings):
    N_bound_0 = inputs["rho"] * inputs["N_total"]
    k0 = ((inputs["N_total"] - N_bound_0) / N_bound_0) * inputs["k_bind"]
    v0 = k0 * inputs["d_step"]
    k_cyl = k0 * (inputs["d_step"] / inputs["l_site"])
    return {"k_ens_walk_0": k0, "v_walk_0": v0, "k_walk_cyl_0": k_cyl}

def symbolic():
    import sympy as sp
    N_total, rho, k_bind, d_step, l_site = sp.symbols("N_total rho k_bind d_step l_site", positive=True)
    N_bound_0 = rho * N_total
    k0 = sp.simplify(((N_total - N_bound_0) / N_bound_0) * k_bind)
    return {
        "k_ens_walk_0": k0,
        "v_walk_0": k0 * d_step,
        "k_walk_cyl_0": k0 * d_step / l_site,
    }

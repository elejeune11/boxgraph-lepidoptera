import math

def run(inputs, settings):
    k_u = (inputs["beta"] * inputs["k_unbind_0"] / inputs["N_bound"]) * math.exp(
        -inputs["F_ext"] / (inputs["N_bound"] * inputs["F_unbind"])
    )
    return {"k_fil_unbind": k_u, "tau_att": 1.0 / k_u}

def symbolic():
    import sympy as sp
    beta, k_unbind_0, N_bound, F_ext, F_unbind = sp.symbols(
        "beta k_unbind_0 N_bound F_ext F_unbind", positive=True
    )
    k_u = (beta * k_unbind_0 / N_bound) * sp.exp(-F_ext / (N_bound * F_unbind))
    return {"k_fil_unbind": k_u, "tau_att": 1 / k_u}

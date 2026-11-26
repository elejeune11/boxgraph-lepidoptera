import math


def run(inputs, settings):
    # With epsilon printed as 0.2 /pN, output dimension is 1/(pN*s) rather than 1/s.
    eps = inputs["epsilon"]
    k0 = inputs["k0_NMIIA_unbind"]
    nb = inputs["N_bound"]
    fext = inputs["F_ext"]
    funb = inputs["F_NMIIA_unbind"]
    return {"k_fil_unbind_printed": (eps * k0 / nb) * math.exp(-fext / (nb * funb))}


def symbolic():
    import sympy as sp
    epsilon, k0_NMIIA_unbind, N_bound, F_ext, F_NMIIA_unbind = sp.symbols(
        "epsilon k0_NMIIA_unbind N_bound F_ext F_NMIIA_unbind", positive=True
    )
    return {
        "k_fil_unbind_printed": (epsilon * k0_NMIIA_unbind / N_bound)
        * sp.exp(-F_ext / (N_bound * F_NMIIA_unbind))
    }

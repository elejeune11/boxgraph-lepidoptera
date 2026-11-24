def run(inputs, settings):
    # Supplementary Table 1 expresses bonded energies in k_B T and pair eps_VV = 1.66 in LJ units;
    # thus T* = k_B T / eps_VV = 1 / eps_VV_lj.
    return {"T_star": 1.0 / inputs["eps_VV_lj"]}


def symbolic():
    import sympy as sp

    eps_VV_lj = sp.Symbol("eps_VV_lj", positive=True)
    return {"T_star": 1 / eps_VV_lj}

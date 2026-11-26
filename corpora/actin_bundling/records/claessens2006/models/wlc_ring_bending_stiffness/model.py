# k_B = 1.380649e-23 J/K = 1.380649e-5 pN*um/K
K_B_PN_UM_PER_K = 1.380649e-5


def run(inputs, settings):
    kappa_B = (
        inputs["c_ring"]
        * (inputs["R_c"] ** 3)
        * K_B_PN_UM_PER_K
        * inputs["T"]
        / (inputs["D_B"] ** 2)
    )
    return {"kappa_B": kappa_B}


def symbolic():
    import sympy as sp

    c_ring, R_c, D_B, T = sp.symbols("c_ring R_c D_B T", positive=True)
    return {"kappa_B": c_ring * R_c**3 * K_B_PN_UM_PER_K * T / D_B**2}

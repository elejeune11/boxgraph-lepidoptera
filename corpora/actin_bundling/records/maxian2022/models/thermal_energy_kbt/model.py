# SI exact Boltzmann constant: 1.380649e-23 J/K = 1.380649e-5 pN*um/K
K_B_PN_UM_PER_K = 1.380649e-5


def run(inputs, settings):
    return {"kBT": K_B_PN_UM_PER_K * inputs["T"]}


def symbolic():
    import sympy as sp

    T = sp.Symbol("T", positive=True)
    return {"kBT": sp.Float(K_B_PN_UM_PER_K) * T}

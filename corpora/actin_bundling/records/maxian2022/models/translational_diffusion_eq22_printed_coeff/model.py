# 1 Pa*s = 1 pN*s/um^2, so (pN*um) / ((Pa*s) * um) = um^2/s directly.


def run(inputs, settings):
    D_t = (inputs["kBT"] / (3.0 * inputs["mu"] * inputs["L"])) * inputs["c_tr"]
    return {"D_t": D_t, "six_D_t": 6.0 * D_t}


def symbolic():
    import sympy as sp

    kBT, mu, L, c_tr = sp.symbols("kBT mu L c_tr", positive=True)
    D_t = kBT * c_tr / (3 * mu * L)
    return {"D_t": D_t, "six_D_t": 6 * D_t}

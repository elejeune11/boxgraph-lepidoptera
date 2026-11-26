# 1 Pa*s = 1 pN*s/um^2, so (pN*um) / ((Pa*s) * um) = um^2/s directly.


def run(inputs, settings):
    c_tr = 3.0 * inputs["alpha"] + inputs["beta"]
    D_t = (inputs["kBT"] / (3.0 * inputs["mu"] * inputs["L"])) * c_tr
    return {"D_t": D_t, "six_D_t": 6.0 * D_t}


def symbolic():
    import sympy as sp

    kBT, mu, L, alpha, beta = sp.symbols("kBT mu L alpha beta", positive=True)
    D_t = kBT * (3 * alpha + beta) / (3 * mu * L)
    return {"D_t": D_t, "six_D_t": 6 * D_t}

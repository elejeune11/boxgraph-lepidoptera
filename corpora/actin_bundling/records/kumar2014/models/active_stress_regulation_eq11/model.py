def run(inputs, settings):
    f_0 = inputs["f_0"]
    beta = inputs["beta"]
    A = inputs["A_norm"]
    I = inputs["I_norm"]
    xi = inputs["xi"]
    f = f_0 + (1.0 + beta) * A / (A + xi) + (1.0 - beta) * I / (I + xi)
    A0_fA = (1.0 + beta) * xi / (A + xi) ** 2
    I0_fI = (1.0 - beta) * xi / (I + xi) ** 2
    return {
        "f": f,
        "A0_fA": A0_fA,
        "I0_fI": I0_fI,
        "sum_sens": A0_fA + I0_fI,
    }


def symbolic():
    import sympy as sp

    f_0, beta = sp.symbols("f_0 beta", real=True)
    A_norm, I_norm, xi = sp.symbols("A_norm I_norm xi", positive=True)
    f = f_0 + (1 + beta) * A_norm / (A_norm + xi) + (1 - beta) * I_norm / (I_norm + xi)
    A0_fA = (1 + beta) * xi / (A_norm + xi) ** 2
    I0_fI = (1 - beta) * xi / (I_norm + xi) ** 2
    return {"f": f, "A0_fA": A0_fA, "I0_fI": I0_fI, "sum_sens": A0_fA + I0_fI}

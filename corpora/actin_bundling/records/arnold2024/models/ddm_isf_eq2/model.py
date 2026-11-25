def run(inputs, settings):
    B = inputs["D_0"]
    A = inputs["D_inf"] - B
    F = 1.0 - (inputs["D_t"] - B) / A
    return {"A": A, "B": B, "F": F}


def symbolic():
    import sympy as sp

    D_t, D_0, D_inf = sp.symbols("D_t D_0 D_inf", real=True)
    B = D_0
    A = D_inf - D_0
    F = 1 - (D_t - B) / A
    return {"A": A, "B": B, "F": F}

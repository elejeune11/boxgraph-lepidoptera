def run(inputs, settings):
    tau = float(inputs["tau"])
    T = float(inputs["T"])
    v_0 = float(inputs["v_0"])
    v_1 = float(inputs["v_1"])
    B = float(inputs["B"])
    D0 = T + tau * v_0**2
    K_r = -B / D0
    K1_r = B * tau * v_0 * v_1 / (D0**2)
    return {"K_r": K_r, "K1_r": K1_r, "lambda_r": K1_r}


def symbolic():
    import sympy as sp

    tau, T, v_0, v_1, B = sp.symbols("tau T v_0 v_1 B", real=True)
    D0 = T + tau * v_0**2
    K_r = -B / D0
    K1_r = B * tau * v_0 * v_1 / (D0**2)
    return {"K_r": K_r, "K1_r": K1_r, "lambda_r": K1_r}

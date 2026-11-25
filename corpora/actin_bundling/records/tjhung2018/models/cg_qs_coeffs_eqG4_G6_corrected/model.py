def run(inputs, settings):
    tau = float(inputs["tau"])
    T = float(inputs["T"])
    v_0 = float(inputs["v_0"])
    v_1 = float(inputs["v_1"])
    gamma_0 = float(inputs["gamma_0"])
    gamma_1 = float(inputs["gamma_1"])
    beta_0 = float(inputs["beta_0"])

    D0 = T + tau * v_0**2
    K_qs = -tau * v_0 * gamma_0 / D0
    K1_qs = -(tau / (2.0 * D0)) * (v_0 * gamma_1 + v_1 * gamma_0 - 2.0 * v_0**2 * v_1 * gamma_0 * tau / D0)
    lambda_qs = tau * v_0 * beta_0 / D0 + K1_qs
    return {"K_qs": K_qs, "K1_qs": K1_qs, "lambda_qs": lambda_qs}


def symbolic():
    import sympy as sp

    tau, T, v_0, v_1, gamma_0, gamma_1, beta_0 = sp.symbols(
        "tau T v_0 v_1 gamma_0 gamma_1 beta_0", real=True
    )
    D0 = T + tau * v_0**2
    K_qs = -tau * v_0 * gamma_0 / D0
    K1_qs = -(tau / (2 * D0)) * (v_0 * gamma_1 + v_1 * gamma_0 - 2 * v_0**2 * v_1 * gamma_0 * tau / D0)
    lambda_qs = tau * v_0 * beta_0 / D0 + K1_qs
    return {"K_qs": K_qs, "K1_qs": K1_qs, "lambda_qs": lambda_qs}

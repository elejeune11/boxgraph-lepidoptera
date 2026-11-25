def run(inputs, settings):
    tau = float(inputs["tau"])
    T = float(inputs["T"])
    v_0 = float(inputs["v_0"])
    v_1 = float(inputs["v_1"])
    gamma_0 = float(inputs["gamma_0"])
    B = float(inputs["B"])
    C = float(inputs["C"])
    denom = (T + tau * v_0**2) ** 2
    return {"zeta_cg": (2.0 * tau * v_0 / denom) * (gamma_0 * C - v_1 * B)}


def symbolic():
    import sympy as sp

    tau, T, v_0, v_1, gamma_0, B, C = sp.symbols("tau T v_0 v_1 gamma_0 B C", real=True)
    return {"zeta_cg": (2 * tau * v_0 / (T + tau * v_0**2) ** 2) * (gamma_0 * C - v_1 * B)}

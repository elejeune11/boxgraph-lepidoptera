K_B = 1.380649e-23  # J / K


def run(inputs, settings):
    T = inputs["T"]
    eta = inputs["eta"]
    a_m = inputs["a"] * 1e-6
    gamma_dot = inputs["gamma_dot"]
    D_c_m2s = K_B * T / (16.0 * eta * a_m)
    D_c_um2s = D_c_m2s * 1e12
    Pe_c = gamma_dot * (inputs["a"] ** 2) / D_c_um2s
    return {"D_c": float(D_c_um2s), "Pe_c": float(Pe_c)}


def symbolic():
    import sympy as sp
    T, eta, a, gamma_dot = sp.symbols("T eta a gamma_dot", positive=True)
    D_c = K_B * T * 1e18 / (16 * eta * a)
    Pe_c = gamma_dot * a**2 / D_c
    return {"D_c": D_c, "Pe_c": Pe_c}

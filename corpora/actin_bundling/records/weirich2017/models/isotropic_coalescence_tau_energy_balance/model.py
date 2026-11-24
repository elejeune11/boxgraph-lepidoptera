def run(inputs, settings):
    # eta [Pa*s] / gamma_I [nN/m = 1e-9 N/m] = 1e9 * (eta/gamma_I) [s/m] = 1e3 * (eta/gamma_I) [s/um]
    # V / (r * L_i^2) is dimensionless; times delta_L [um] gives [s]
    inv_vcap_s_per_um = 1e3 * inputs["eta"] / inputs["gamma_I"]
    geom = inputs["V"] / (inputs["r"] * (inputs["L_i"] ** 2))
    return {"tau": inv_vcap_s_per_um * geom * inputs["delta_L"]}


def symbolic():
    import sympy as sp

    eta, gamma_I, V, r, L_i, delta_L = sp.symbols("eta gamma_I V r L_i delta_L", positive=True)
    return {"tau": sp.Integer(1000) * (eta / gamma_I) * (V / (r * L_i**2)) * delta_L}

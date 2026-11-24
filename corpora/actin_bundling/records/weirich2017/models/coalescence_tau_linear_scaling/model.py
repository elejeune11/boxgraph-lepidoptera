def run(inputs, settings):
    # eta [Pa*s] / gamma_eff [nN/m] = 1e3 * eta / gamma_eff [s/um]
    slope = 1e3 * inputs["eta"] / inputs["gamma_eff"]
    return {"slope": slope, "tau": slope * inputs["delta_L"]}


def symbolic():
    import sympy as sp

    eta, gamma_eff, delta_L = sp.symbols("eta gamma_eff delta_L", positive=True)
    slope = sp.Integer(1000) * eta / gamma_eff
    return {"slope": slope, "tau": slope * delta_L}

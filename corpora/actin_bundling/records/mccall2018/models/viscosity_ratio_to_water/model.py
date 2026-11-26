def run(inputs, settings):
    # eta in Pa*s, eta_water in mPa*s (1 Pa*s = 1000 mPa*s)
    return {"eta_ratio": (inputs["eta"] * 1000.0) / inputs["eta_water"]}


def symbolic():
    import sympy as sp

    eta, eta_water = sp.symbols("eta eta_water", positive=True)
    return {"eta_ratio": 1000 * eta / eta_water}

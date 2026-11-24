def run(inputs, settings):
    # sigma in nm -> 1e-9 m; m in kg; eps_VV in pN*nm -> 1e-21 J (= kg*m^2/s^2)
    # sigma * m / eps_VV has units m * kg / (kg*m^2/s^2) = s^2 / m (NOT time!)
    sigma_m = inputs["sigma"] * 1e-9
    eps_J = inputs["eps_VV"] * 1e-21
    return {"tau_MD_printed": sigma_m * (inputs["m"] / eps_J)}


def symbolic():
    import sympy as sp

    sigma, m, eps_VV = sp.symbols("sigma m eps_VV", positive=True)
    return {"tau_MD_printed": (sigma * sp.Float("1e-9")) * (m / (eps_VV * sp.Float("1e-21")))}

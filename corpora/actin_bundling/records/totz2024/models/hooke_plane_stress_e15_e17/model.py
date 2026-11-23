def run(inputs, settings):
    Y = inputs["Y"]
    nu = inputs["nu"]
    exx = inputs["eps_xx"]
    eyy = inputs["eps_yy"]
    exy = inputs["eps_xy"]
    pref = Y / (1.0 - nu ** 2)
    return {
        "sigma_xx": pref * (exx + nu * eyy),
        "sigma_yy": pref * (eyy + nu * exx),
        "sigma_xy": (Y / (1.0 + nu)) * exy,
        "eps_zz": -(nu / (1.0 - nu)) * (exx + eyy),
    }


def symbolic():
    import sympy as sp
    Y, nu, eps_xx, eps_yy, eps_xy = sp.symbols("Y nu eps_xx eps_yy eps_xy", positive=True)
    return {
        "sigma_xx": Y / (1 - nu**2) * (eps_xx + nu * eps_yy),
        "sigma_yy": Y / (1 - nu**2) * (eps_yy + nu * eps_xx),
        "sigma_xy": Y / (1 + nu) * eps_xy,
        "eps_zz": -nu / (1 - nu) * (eps_xx + eps_yy),
    }

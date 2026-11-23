def run(inputs, settings):
    Y = inputs["Y"]
    h = inputs["h"]
    nu = inputs["nu"]
    D = Y * (h ** 3) / (12.0 * (1.0 - nu ** 2))
    return {"D": D}


def symbolic():
    import sympy as sp
    Y, h, nu = sp.symbols("Y h nu", positive=True)
    return {"D": Y * h**3 / (12 * (1 - nu**2))}

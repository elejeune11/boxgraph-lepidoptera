def run(inputs, settings):
    R = float(inputs["R"])
    x = float(inputs["x"])
    return {"S": R * R / x, "V": (R**3) / (x * x)}


def symbolic():
    import sympy as sp

    R, x = sp.symbols("R x", positive=True)
    return {"S": R**2 / x, "V": R**3 / x**2}

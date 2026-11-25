def run(inputs, settings):
    return {"ratio": inputs["eta_droplet"] / inputs["eta_water"]}


def symbolic():
    import sympy as sp

    eta_droplet, eta_water = sp.symbols("eta_droplet eta_water", positive=True)
    return {"ratio": eta_droplet / eta_water}

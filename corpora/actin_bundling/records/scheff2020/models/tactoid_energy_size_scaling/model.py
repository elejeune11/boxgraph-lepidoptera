def run(inputs, settings):
    R = float(inputs["R"])
    return {"F_s": float(inputs["gamma"]) * R * R, "F_el": float(inputs["K"]) * R}


def symbolic():
    import sympy as sp

    R, K, gamma = sp.symbols("R K gamma", positive=True)
    return {"F_s": gamma * R**2, "F_el": K * R}

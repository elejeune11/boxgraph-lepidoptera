def run(inputs, settings):
    # 2a = L and a*b = R^2 => a = L/2, b = 2*R^2/L => a/b = L^2 / (4 * R^2)
    return {"ar_max": (inputs["L"] ** 2) / (4.0 * (inputs["R"] ** 2))}


def symbolic():
    import sympy as sp

    L, R = sp.symbols("L R", positive=True)
    return {"ar_max": L ** 2 / (4 * R ** 2)}

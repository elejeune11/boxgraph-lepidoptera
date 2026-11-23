def run(inputs, settings):
    # eta_3D [Pa*s] * h [nm -> m = 1e-9] -> eta [Pa*s*m]
    return {"eta": inputs["eta_3D"] * (inputs["h"] * 1e-9)}


def symbolic():
    import sympy as sp

    eta_3D, h = sp.symbols("eta_3D h", positive=True)
    return {"eta": eta_3D * h * sp.Rational(1, 10**9)}

import math


def run(inputs, settings):
    theta = float(inputs["theta"])
    L = float(inputs["L"])          # nm
    D = float(inputs["D"])          # nm
    L_0 = float(inputs["L_0"])      # nm
    L_p = float(inputs["L_p"])      # um
    kBT = float(inputs["kBT"])      # pN*nm

    if abs(theta) < 1e-12:
        x = 0.0
        y = 0.0
    else:
        x = L * (1.0 - math.cos(theta)) / theta
        y = L * (1.0 - math.sin(theta) / theta)

    d = D - 2.0 * x
    L_prot = L - L_0 - y
    # L_p is in um -> convert to nm via * 1e3; curvature kappa_c = theta / L
    E_fil = kBT * (L_p * 1e3 / (2.0 * L)) * (theta ** 2)
    return {
        "x": x,
        "y": y,
        "d": d,
        "L_prot": L_prot,
        "E_fil": E_fil,
    }


def symbolic():
    import sympy as sp

    theta, L, D, L_0, L_p, kBT = sp.symbols("theta L D L_0 L_p kBT", positive=True)
    x = L * (1 - sp.cos(theta)) / theta
    y = L * (1 - sp.sin(theta) / theta)
    d = D - 2 * x
    L_prot = L - L_0 - y
    E_fil = kBT * (L_p * 1000 / (2 * L)) * theta**2
    return {
        "x": x,
        "y": y,
        "d": d,
        "L_prot": L_prot,
        "E_fil": E_fil,
    }

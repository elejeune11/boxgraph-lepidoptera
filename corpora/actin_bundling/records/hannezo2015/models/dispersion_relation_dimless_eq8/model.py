def run(inputs, settings):
    q = inputs["q_hat"]
    chi = inputs["chi"]
    phi = inputs["phi"]
    q2 = q * q
    lam_hat = -phi - q2 + chi * q2 / (1.0 + q2)
    det_L = q2 * q2 + (1.0 + phi - chi) * q2 + phi
    tr_L = 2.0 * q2 + 1.0 + phi - chi
    return {"lam_hat": lam_hat, "det_L": det_L, "tr_L": tr_L}


def symbolic():
    import sympy as sp

    q_hat, chi, phi = sp.symbols("q_hat chi phi", real=True)
    q2 = q_hat**2
    return {
        "lam_hat": -phi - q2 + chi * q2 / (1 + q2),
        "det_L": q2**2 + (1 + phi - chi) * q2 + phi,
        "tr_L": 2 * q2 + 1 + phi - chi,
    }

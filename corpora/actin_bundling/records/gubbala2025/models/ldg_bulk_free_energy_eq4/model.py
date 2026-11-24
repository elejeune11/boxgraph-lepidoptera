def run(inputs, settings):
    beta = float(inputs["beta"])
    phi = float(inputs["phi"])
    Q_xx = float(inputs["Q_xx"])
    Q_xy = float(inputs["Q_xy"])
    q2 = Q_xx * Q_xx + Q_xy * Q_xy
    # In 2D for traceless symmetric Q: Tr(Q^2) = 2*q2, Tr(Q^4) = 2*q2^2
    tr_Q2 = 2.0 * q2
    tr_Q4 = 2.0 * q2 * q2
    f_bulk = 0.25 * beta * tr_Q4 - 0.5 * (1.0 + phi) * tr_Q2
    return {"f_bulk": f_bulk, "tr_Q2": tr_Q2}


def symbolic():
    import sympy as sp

    beta, phi, Q_xx, Q_xy = sp.symbols("beta phi Q_xx Q_xy", real=True)
    q2 = Q_xx**2 + Q_xy**2
    tr_Q2 = 2 * q2
    tr_Q4 = 2 * q2**2
    return {
        "f_bulk": beta * tr_Q4 / 4 - (1 + phi) * tr_Q2 / 2,
        "tr_Q2": tr_Q2,
    }

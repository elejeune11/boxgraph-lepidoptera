import cmath

def run(inputs, settings):
    nu = inputs["nu"]
    x = nu * nu
    kd = inputs["k_d"]
    L = inputs["L"]
    c0 = inputs["c_0"]
    b = inputs["b"]
    er = inputs["eta_rot"]
    beta = inputs["beta"]
    lo = inputs["lambda_o"]
    kap = inputs["kappa"]
    lam = inputs["lam"]
    q0 = inputs["q_0"]

    G = c0 if c0 >= 0.0 else 8.0 * b * q0 * q0
    f_vis = 4.0 + 1.0 / x
    g_nem = 2.0 * G + 2.0 * L * x
    A = 2.0 * er * f_vis - beta * beta
    B = (
        f_vis * g_nem
        + 2.0 * er * (f_vis * (x + kd) - lam * (1.0 + kap * q0))
        - beta * (beta * (x + kd) + 2.0 * lo * q0 + lam * kap)
    )
    C = (f_vis * (x + kd) - lam * (1.0 + kap * q0)) * g_nem - lam * kap * (
        beta * (x + kd) + 2.0 * lo * q0
    )
    disc = cmath.sqrt(B * B - 4.0 * A * C)
    ap = ((-B + disc) / (2.0 * A)).real
    am = ((-B - disc) / (2.0 * A)).real
    return {"alpha_plus": max(ap, am), "alpha_minus": min(ap, am)}

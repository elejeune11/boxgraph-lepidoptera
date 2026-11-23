def run(inputs, settings):
    rho_a = inputs["rho_a"]
    rho_s = inputs["rho_s"]
    if rho_a <= 0.0 or rho_s <= 0.0:
        raise ValueError("Reaction rate constants rho_a and rho_s must be positive")
    # Solve R_a = rho_a * (a^2 * s - a) = 0 and R_s = rho_s * (1 - a^2 * s) = 0
    # by 2D Newton iteration from an off-equilibrium positive guess.
    a, s = 0.5, 2.0
    for _ in range(50):
        ra = rho_a * (a * a * s - a)
        rs = rho_s * (1.0 - a * a * s)
        if abs(ra) < 1e-14 and abs(rs) < 1e-14:
            break
        j11 = rho_a * (2.0 * a * s - 1.0)
        j12 = rho_a * (a * a)
        j21 = rho_s * (-2.0 * a * s)
        j22 = rho_s * (-a * a)
        det = j11 * j22 - j12 * j21
        da = -(j22 * ra - j12 * rs) / det
        ds = -(-j21 * ra + j11 * rs) / det
        a += da
        s += ds
    return {"a0": float(a), "s0": float(s)}


def symbolic():
    import sympy as sp

    rho_a, rho_s = sp.symbols("rho_a rho_s", positive=True)
    # Unique positive fixed point of rho_a*(a^2*s - a) = 0 and rho_s*(1 - a^2*s) = 0
    return {"a0": sp.Integer(1) + 0 * rho_a, "s0": sp.Integer(1) + 0 * rho_s}

import math


def run(inputs, settings):
    # Axisymmetric 3D distribution with <u_x^2> = cos^2(theta), <u_y^2> = <u_z^2> = sin^2(theta)/2.
    # Printed Eq. (S18) sums i=0..N (N+1 terms) while dividing by N, multiplying the average by (N+1)/N.
    N = inputs["N"]
    c2 = math.cos(inputs["theta"]) ** 2
    s2 = 1.0 - c2
    scale = (N + 1.0) / N
    q_xx = 1.5 * (scale * c2 - 1.0 / 3.0)
    q_perp = 1.5 * (scale * (0.5 * s2) - 1.0 / 3.0)
    return {"S": max(q_xx, q_perp)}


def symbolic():
    import sympy as sp
    N, theta = sp.symbols("N theta", positive=True)
    # For aligned/nematic branch (cos^2(theta) >= 1/3), max eigenvalue is Q_xx:
    return {"S": sp.Rational(3, 2) * (((N + 1) / N) * sp.cos(theta) ** 2 - sp.Rational(1, 3))}

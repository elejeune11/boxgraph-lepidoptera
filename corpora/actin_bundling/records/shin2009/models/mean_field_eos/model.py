import math
from scipy.optimize import brentq


def grand_potential(rho: float, z: float, u: float) -> float:
    return (
        -0.5 * u * rho * rho
        + rho * math.log(rho)
        + (1.0 - rho) * math.log(1.0 - rho)
        - rho * math.log(z)
    )


def solve_equilibrium_rho(z: float, u: float, tol: float) -> float:
    """Find all roots of f(rho) = rho - 1/(1 + exp(-u*rho)/z) = 0 on (0, 1)
    and select the global minimum of the mean-field grand potential."""
    def f(r):
        return r - 1.0 / (1.0 + math.exp(-u * r) / z)

    eps = 1.0e-12
    # Partition (eps, 1-eps) around spinodal points where 1 - u*r*(1-r) = 0
    cuts = [eps]
    if u > 4.0:
        disc = math.sqrt(1.0 - 4.0 / u)
        cuts.extend([0.5 * (1.0 - disc), 0.5 * (1.0 + disc)])
    cuts.append(1.0 - eps)

    roots = []
    for a, b in zip(cuts[:-1], cuts[1:]):
        fa, fb = f(a), f(b)
        if fa == 0.0:
            roots.append(a)
        elif fa * fb < 0.0:
            roots.append(brentq(f, a, b, xtol=tol, rtol=tol))
        elif fb == 0.0:
            roots.append(b)

    if not roots:
        # Fallback fine bracket scan
        n_grid = 200
        grid = [eps + (1.0 - 2.0 * eps) * i / n_grid for i in range(n_grid + 1)]
        for a, b in zip(grid[:-1], grid[1:]):
            if f(a) * f(b) <= 0.0:
                roots.append(brentq(f, a, b, xtol=tol, rtol=tol))

    # Choose thermodynamically stable root minimizing grand potential
    best_rho = min(roots, key=lambda r: grand_potential(r, z, u))
    return float(best_rho)


def run(inputs, settings):
    z = float(inputs["z"])
    u = float(inputs["u"])
    tol = float(settings.get("tol", 1.0e-10))

    rho = solve_equilibrium_rho(z, u, tol)
    denom = z * (1.0 - u * rho * (1.0 - rho))
    chi_rho = (rho * (1.0 - rho) / denom) if abs(denom) > 1.0e-30 else float("inf")

    z_coex = math.exp(-0.5 * u)
    if u > 4.0:
        # Symmetric coexistence: ln(r/(1-r)) = u*(r - 1/2)
        rho_hi = brentq(
            lambda r: math.log(r / (1.0 - r)) - u * (r - 0.5),
            0.5 + 1.0e-9,
            1.0 - 1.0e-12,
            xtol=tol,
            rtol=tol,
        )
        delta_rho_coex = 2.0 * rho_hi - 1.0
    else:
        delta_rho_coex = 0.0

    return {
        "rho": rho,
        "chi_rho": chi_rho,
        "z_coex": z_coex,
        "delta_rho_coex": delta_rho_coex,
    }

"""Einstein-Stokes scaling of diffusion coefficient for globular proteins of equal density: D ~ R_H^-1 ~ M^(-1/3)."""


def run(inputs, settings):
    d_ref = float(inputs["D_ref"])
    m_ref = float(inputs["M_ref"])
    m = float(inputs["M"])
    return {"D": d_ref * (m_ref / m) ** (1.0 / 3.0)}


def symbolic():
    import sympy as sp

    D_ref, M_ref, M = sp.symbols("D_ref M_ref M", positive=True)
    return {"D": D_ref * (M_ref / M) ** sp.Rational(1, 3)}

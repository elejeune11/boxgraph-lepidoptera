import math

def run(inputs, settings):
    # 1 pN = 1e-12 N = 1 Pa * um^2, so L [pN] / c0 [Pa] has units um^2.
    abs_c0 = abs(2.0 * inputs["a"] - inputs["rho_lam_o"])
    ell_q = math.sqrt(inputs["L"] / abs_c0)
    t_q = inputs["eta_rot"] / abs_c0
    return {"ell_q": ell_q, "t_q": t_q}

def symbolic():
    import sympy as sp
    L, a, rho_lam_o, eta_rot = sp.symbols("L a rho_lam_o eta_rot", positive=True)
    abs_c0 = sp.Abs(2 * a - rho_lam_o)
    return {"ell_q": sp.sqrt(L / abs_c0), "t_q": eta_rot / abs_c0}

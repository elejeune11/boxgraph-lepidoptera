def run(inputs, settings):
    kappa = inputs["kappa"]
    kBT = inputs["kBT"]
    Lambda = inputs["Lambda"]
    H0 = inputs["H0"]
    ns = inputs["ns"]
    J = inputs["J"]
    n = inputs["n"]
    grad_n = inputs["grad_n"]
    grad3_n = inputs["grad3_n"]
    grad_H = inputs["grad_H"]

    J_diff = -Lambda * kBT * grad_n
    J_disp = -(kappa * Lambda * (H0 ** 2) / ns) * n * grad_n
    J_agg = (2.0 * J * Lambda / ns) * n * (grad_n + grad3_n / ns)
    J_curv = (kappa * Lambda * H0 / ns) * n * grad_H
    return {"J_tot": J_diff + J_disp + J_agg + J_curv}


def symbolic():
    import sympy as sp

    kappa, kBT, Lambda, H0, ns, J, n, grad_n, grad3_n, grad_H = sp.symbols(
        "kappa kBT Lambda H0 ns J n grad_n grad3_n grad_H", real=True
    )
    J_diff = -Lambda * kBT * grad_n
    J_disp = -(kappa * Lambda * H0 ** 2 / ns) * n * grad_n
    J_agg = (2 * J * Lambda / ns) * n * (grad_n + grad3_n / ns)
    J_curv = (kappa * Lambda * H0 / ns) * n * grad_H
    return {"J_tot": J_diff + J_disp + J_agg + J_curv}

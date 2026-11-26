"""2D Hermitian beam FE model with randomly aligned fiber fractures (Bathe et al. 2008, p. 5 & SM p. 11)."""

import math
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


def run(inputs, settings):
    M = int(round(float(inputs["M"])))
    alpha = float(inputs["alpha"])
    phi = float(inputs["phi"])
    kappa_f_sim = float(inputs["kappa_f_sim"])
    n_el = int(settings.get("n_el", 80))
    seed = int(settings.get("seed", 0))

    rng = np.random.default_rng(seed)
    n_fibers = 2 * M
    L = 1.0
    df = 1e-3
    Ef = 1.0
    EA = Ef * (df ** 2)
    EI = Ef * (df ** 4) / 12.0
    kx_over_delta = alpha * EA / (L ** 2)

    n_el_eff = max(n_el, int(math.ceil(4.0 * phi)))
    if n_el_eff % 2 == 1:
        n_el_eff += 1
    n_nodes = n_el_eff + 1
    h = L / n_el_eff

    cuts = []
    for _ in range(n_fibers):
        if phi <= 1.0:
            cuts.append(set())
        else:
            phase = rng.uniform(0.0, L / phi)
            xs = [phase + m * (L / phi) for m in range(int(math.ceil(phi)) + 1)]
            cut_nodes = {int(round(x / h)) for x in xs if 0 < int(round(x / h)) < n_el_eff}
            cuts.append(cut_nodes)

    u_left = np.zeros((n_fibers, n_el_eff), dtype=int)
    u_right = np.zeros((n_fibers, n_el_eff), dtype=int)
    th_left = np.zeros((n_fibers, n_el_eff), dtype=int)
    th_right = np.zeros((n_fibers, n_el_eff), dtype=int)

    next_dof = n_nodes
    for k in range(n_fibers):
        cur_u = next_dof
        next_dof += 1
        cur_th = next_dof
        next_dof += 1
        for e in range(n_el_eff):
            u_left[k, e] = cur_u
            th_left[k, e] = cur_th
            cur_u = next_dof
            next_dof += 1
            cur_th = next_dof
            next_dof += 1
            u_right[k, e] = cur_u
            th_right[k, e] = cur_th
            if (e + 1) in cuts[k]:
                cur_u = next_dof
                next_dof += 1
                cur_th = next_dof
                next_dof += 1

    n_dof = next_dof
    rows, cols, vals = [], [], []

    def add_mat(dofs, Ke):
        for a, da in enumerate(dofs):
            for b, db in enumerate(dofs):
                rows.append(da)
                cols.append(db)
                vals.append(Ke[a, b])

    ka = (EA / h) * np.array([[1.0, -1.0], [-1.0, 1.0]])
    kb = (EI / (h ** 3)) * np.array([
        [12.0, 6.0 * h, -12.0, 6.0 * h],
        [6.0 * h, 4.0 * h * h, -6.0 * h, 2.0 * h * h],
        [-12.0, -6.0 * h, 12.0, -6.0 * h],
        [6.0 * h, 2.0 * h * h, -6.0 * h, 4.0 * h * h],
    ])
    for k in range(n_fibers):
        for e in range(n_el_eff):
            add_mat([u_left[k, e], u_right[k, e]], ka)
            add_mat([e, th_left[k, e], e + 1, th_right[k, e]], kb)

    c_th = 0.5 * df
    v = np.array([1.0, -1.0, c_th, c_th])
    K_cl_unit = np.outer(v, v)
    kx_half = 0.5 * kx_over_delta * h
    for k in range(1, n_fibers):
        for e in range(n_el_eff):
            add_mat([u_left[k, e], u_left[k - 1, e], th_left[k, e], th_left[k - 1, e]], kx_half * K_cl_unit)
            add_mat([u_right[k, e], u_right[k - 1, e], th_right[k, e], th_right[k - 1, e]], kx_half * K_cl_unit)

    K = sp.coo_matrix((vals, (rows, cols)), shape=(n_dof, n_dof)).tocsc()
    f = np.zeros(n_dof)
    mid = n_el_eff // 2
    f[mid] = 1.0
    fixed = {0, n_el_eff, int(u_left[0, mid])}
    free = np.array([i for i in range(n_dof) if i not in fixed])
    u = np.zeros(n_dof)
    u[free] = spla.spsolve(K[free][:, free], f[free])
    w_mid = float(u[mid])
    kappa_2d = (1.0 * (L ** 3)) / (48.0 * w_mid)
    # 3D normalized kappa_B* = 2M * kappa_2d / EI = (n_fibers) * kappa_2d / EI
    # In Fig. S1, unnormalized 2D stiffness kappa_B_sim = kappa_f_sim * (kappa_2d / EI)
    kappa_2d_over_EI = kappa_2d / EI
    return {
        "kappa_B_star": float(n_fibers * kappa_2d_over_EI),
        "kappa_B_sim": float(kappa_f_sim * kappa_2d_over_EI),
    }

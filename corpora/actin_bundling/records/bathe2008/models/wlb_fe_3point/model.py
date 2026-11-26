"""2D Hermitian beam + discrete cross-link shear Finite Element model (Bathe et al. 2008, p. 3).

Implements the exact finite element formulation described on p. 3 col 1-2:
- 2M = sqrt(N) planar Euler-Bernoulli fibers discretized into n_el two-node Hermitian beam elements {u_i, w_i, theta_i}
- Adjacent fibers constrained to equal transverse deflection w_i
- Nearest-neighbor cross-link shear energy E = (k_x/2)[(u^(k)-u^(k-1)) + ((d_f+t)/2)*(theta^(k)+theta^(k-1))]^2
- Three-point bending midpoint unit load P = 1 with pinned (a = 48) or clamped (a = 192) ends
- 3D bundle stiffness kappa_B = 2M * kappa_{B(2D)}, normalized by kappa_f = E_f * I_f.
"""

import math
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


def run(inputs, settings):
    N = float(inputs["N"])
    alpha = float(inputs["alpha"])
    chi2 = float(inputs["chi2"])
    bc = settings.get("bc", "pinned")
    n_el = int(settings.get("n_el", 80))
    if n_el % 2 == 1:
        n_el += 1

    # Number of planar layers 2M = sqrt(N) (integer in the discrete FE model; if N is not a perfect square, round)
    n_fibers = max(1, int(round(math.sqrt(N))))
    if n_fibers == 1:
        return {"kappa_B_star": 1.0}

    n_nodes = n_el + 1
    L = 1.0
    h = L / n_el
    df = 1e-3
    # Match chi^2 = A_f * (d_f + t)^2 / (12 * I_f) using square fiber A_f = df^2, I_f = df^4/12 => (df+t) = df * sqrt(chi^2)
    d_eff = df * math.sqrt(chi2)
    Ef = 1.0
    EA = Ef * (df ** 2)
    EI = Ef * (df ** 4) / 12.0
    kx_over_delta = alpha * EA / (L ** 2)

    n_dof = n_nodes + 2 * n_fibers * n_nodes

    def idx_w(j):
        return j

    def idx_u(k, j):
        return n_nodes + k * 2 * n_nodes + j

    def idx_th(k, j):
        return n_nodes + k * 2 * n_nodes + n_nodes + j

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
        for e in range(n_el):
            add_mat([idx_u(k, e), idx_u(k, e + 1)], ka)
            add_mat([idx_w(e), idx_th(k, e), idx_w(e + 1), idx_th(k, e + 1)], kb)

    c_th = 0.5 * d_eff
    v = np.array([1.0, -1.0, c_th, c_th])
    K_cl_unit = np.outer(v, v)
    for j in range(n_nodes):
        wt = 0.5 if (j == 0 or j == n_el) else 1.0
        kx_node = kx_over_delta * h * wt
        for k in range(1, n_fibers):
            add_mat([idx_u(k, j), idx_u(k - 1, j), idx_th(k, j), idx_th(k - 1, j)], kx_node * K_cl_unit)

    K = sp.coo_matrix((vals, (rows, cols)), shape=(n_dof, n_dof)).tocsc()
    f = np.zeros(n_dof)
    mid = n_el // 2
    f[idx_w(mid)] = 1.0

    fixed = {idx_w(0), idx_w(n_el)}
    for k in range(n_fibers):
        fixed.add(idx_u(k, mid))
    if bc == "clamped":
        for k in range(n_fibers):
            fixed.add(idx_th(k, 0))
            fixed.add(idx_th(k, n_el))
            fixed.add(idx_u(k, 0))
            fixed.add(idx_u(k, n_el))

    free = np.array([i for i in range(n_dof) if i not in fixed])
    u = np.zeros(n_dof)
    u[free] = spla.spsolve(K[free][:, free], f[free])
    w_mid = float(u[idx_w(mid)])
    a_fac = 48.0 if bc == "pinned" else 192.0
    kappa_2d = (1.0 * (L ** 3)) / (a_fac * w_mid)
    # Scale by (N / n_fibers) rather than n_fibers if N is non-square (exact 2M when N=(2M)^2)
    kappa_B_star = (N / n_fibers) * kappa_2d / EI
    return {"kappa_B_star": float(kappa_B_star)}

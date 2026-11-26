"""Axisymmetric P1-FEM solver for the finite-pool cone kinetic PDEs Eqs. (6)-(16) as printed.

Implements printed boundary conditions Eqs. (14)-(15) verbatim (without +k_2*V* and +k_4*A*).
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import bmat, coo_matrix, diags

NA_PER_UM3_PER_UM = 602.214076


def _solve_cone(inputs, settings, conserved_bc, surf_axisym):
    G0 = float(inputs["G_0"])
    Arp0 = float(inputs["Arp_0"])
    CP0 = float(inputs["CP_0"]) * 1e-3  # nM -> uM
    V0 = float(inputs["V_0"])
    B0_norm = float(inputs["B_0_norm"])
    sigma = float(inputs["sigma"]) * 1e-3  # nm -> um
    R_cone = float(inputs["r_cone"])
    h = float(inputs["h_cone"])
    kon = float(inputs["k_on"])
    koff = float(inputs["k_off"])
    k1 = float(inputs["k_1"])
    k2 = float(inputs["k_2"])
    k3 = float(inputs["k_3"]) / float(inputs["k_3_reduction"])
    k4 = float(inputs["k_4"])
    kCap = float(inputs["k_Cap"])
    kArp = float(inputs["k_Arp"])
    kre = float(inputs["k_re"])
    kDeg = float(inputs["k_Deg"])
    DG = float(inputs["D_G"])
    DArp = float(inputs["D_Arp"])
    DCP = float(inputs["D_CP"])
    DV = float(inputs["D_V"])
    DVs = float(inputs["D_Vstar"])
    DA = float(inputs["D_A"])
    DB = float(inputs["D_B"])

    mesh_n = max(12, int(round(float(settings.get("mesh_n", 72)))))
    Nr = mesh_n
    Nz = max(8, mesh_n // 3)
    r_stretch = float(settings.get("r_stretch", 3.0))
    z_stretch = float(settings.get("z_stretch", 2.0))
    t_end = float(settings.get("t_end", 40.0))
    rtol = float(settings.get("rtol", 1e-5))
    atol = float(settings.get("atol", 1e-7))

    n_in = max(6, Nr // 4)
    n_out = Nr - n_in
    xi_split = min(0.25 / R_cone, 0.25)
    s_in = np.linspace(0.0, 1.0, n_in + 1)
    xi_in = xi_split * (np.sinh(r_stretch * s_in) / np.sinh(r_stretch) if r_stretch > 0 else s_in)
    xi_out = np.linspace(xi_split, 1.0, n_out + 1)[1:]
    xi = np.concatenate([xi_in, xi_out])
    s_z = np.linspace(0.0, 1.0, Nz + 1)
    zeta = np.sinh(z_stretch * s_z) / np.sinh(z_stretch) if z_stretch > 0 else s_z

    nodes_r, nodes_z = [], []
    node_id = np.full((Nr + 1, Nz + 1), -1, dtype=int)
    idx = 0
    for j in range(Nz + 1):
        z = h * zeta[j]
        if j == Nz:
            nodes_r.append(0.0)
            nodes_z.append(h)
            node_id[:, j] = idx
            idx += 1
        else:
            scale = 1.0 - zeta[j]
            for i in range(Nr + 1):
                nodes_r.append(R_cone * xi[i] * scale)
                nodes_z.append(z)
                node_id[i, j] = idx
                idx += 1
    nodes_r = np.array(nodes_r)
    nodes_z = np.array(nodes_z)
    Nb = len(nodes_r)

    tris = []
    for j in range(Nz):
        for i in range(Nr):
            n00, n10, n01, n11 = node_id[i, j], node_id[i + 1, j], node_id[i, j + 1], node_id[i + 1, j + 1]
            if j == Nz - 1:
                tris.append((n00, n10, n01))
            else:
                tris.append((n00, n10, n11))
                tris.append((n00, n11, n01))
    tris = np.array(tris, dtype=int)

    r0, z0 = nodes_r[tris[:, 0]], nodes_z[tris[:, 0]]
    r1, z1 = nodes_r[tris[:, 1]], nodes_z[tris[:, 1]]
    r2, z2 = nodes_r[tris[:, 2]], nodes_z[tris[:, 2]]
    area = 0.5 * np.abs((r1 - r0) * (z2 - z0) - (r2 - r0) * (z1 - z0))
    rc = (r0 + r1 + r2) / 3.0
    b = np.stack([z1 - z2, z2 - z0, z0 - z1], axis=1)
    c = np.stack([r2 - r1, r0 - r2, r1 - r0], axis=1)

    M_bulk = np.zeros(Nb)
    M_e = (2.0 * np.pi * area[:, None] / 12.0) * np.stack(
        [2.0 * r0 + r1 + r2, r0 + 2.0 * r1 + r2, r0 + r1 + 2.0 * r2], axis=1
    )
    for a in range(3):
        np.add.at(M_bulk, tris[:, a], M_e[:, a])

    rows, cols, vals = [], [], []
    fac = 2.0 * np.pi * rc / (4.0 * area)
    for a in range(3):
        for bb in range(3):
            rows.append(tris[:, a])
            cols.append(tris[:, bb])
            vals.append(fac * (b[:, a] * b[:, bb] + c[:, a] * c[:, bb]))
    K_bulk = coo_matrix(
        (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(Nb, Nb)
    ).tocsc()

    Nm = Nr + 1
    r_m = nodes_r[:Nm]
    dr_m = np.diff(r_m)
    rm_mid = 0.5 * (r_m[:-1] + r_m[1:])
    i0, i1 = np.arange(Nr), np.arange(1, Nm)

    M_surf_ax = np.zeros(Nm)
    np.add.at(M_surf_ax, i0, 2.0 * np.pi * dr_m / 6.0 * (2.0 * r_m[:-1] + r_m[1:]))
    np.add.at(M_surf_ax, i1, 2.0 * np.pi * dr_m / 6.0 * (r_m[:-1] + 2.0 * r_m[1:]))

    if surf_axisym:
        M_s = M_surf_ax
        sf = 2.0 * np.pi * rm_mid / dr_m
    else:
        M_s = np.zeros(Nm)
        np.add.at(M_s, i0, 0.5 * dr_m)
        np.add.at(M_s, i1, 0.5 * dr_m)
        sf = 1.0 / dr_m

    K_s = coo_matrix(
        (np.concatenate([sf, -sf, -sf, sf]),
         (np.concatenate([i0, i0, i1, i1]), np.concatenate([i0, i1, i0, i1]))),
        shape=(Nm, Nm),
    ).tocsc()
    L_s = K_s.multiply(1.0 / M_s[:, None]).tocsc()
    L_b = K_bulk.multiply(1.0 / M_bulk[:, None]).tocsc()
    s2b = (M_surf_ax / M_bulk[:Nm]) / NA_PER_UM3_PER_UM

    L_s_V = -DV * L_s
    L_s_Vs = -DVs * L_s
    L_s_A = -DA * L_s
    L_b_G = -DG * L_b
    L_b_A = -DArp * L_b
    L_b_C = -DCP * L_b
    cb = 1.0 if conserved_bc else 0.0
    N_s = 5 * Nm
    N_tot = N_s + 3 * Nb

    def rhs(t, y):
        V, Vs, As, B, Vp = y[0:Nm], y[Nm:2 * Nm], y[2 * Nm:3 * Nm], y[3 * Nm:4 * Nm], y[4 * Nm:5 * Nm]
        G = y[N_s:N_s + Nb]
        Arp = y[N_s + Nb:N_s + 2 * Nb]
        CP = y[N_s + 2 * Nb:N_tot]
        Gm, Arpm, CPm = G[:Nm], Arp[:Nm], CP[:Nm]

        dV = L_s_V @ V - k1 * Gm * V + k2 * Vs + kre * Vp
        dVs = L_s_Vs @ Vs + k1 * Gm * V - k2 * Vs - k3 * Arpm * Vs + k4 * As
        dAs = L_s_A @ As + k3 * Arpm * Vs - k4 * As - kArp * As * B
        dB = -DB * (Gm / G0) * (L_s @ B) + kArp * As * B - kCap * CPm * B - kDeg * B
        dVp = kArp * As * B - kre * Vp

        dG = L_b_G @ G
        dArp = L_b_A @ Arp
        dCP = L_b_C @ CP
        dG[:Nm] += s2b * (-k1 * Gm * V + cb * k2 * Vs - kon * Gm * B + koff * B)
        dArp[:Nm] += s2b * (-k3 * Arpm * Vs + cb * k4 * As)
        dCP[:Nm] += s2b * (-kCap * CPm * B)
        return np.concatenate([dV, dVs, dAs, dB, dVp, dG, dArp, dCP])

    def jac(t, y):
        V, Vs, As, B, Vp = y[0:Nm], y[Nm:2 * Nm], y[2 * Nm:3 * Nm], y[3 * Nm:4 * Nm], y[4 * Nm:5 * Nm]
        G = y[N_s:N_s + Nb]
        Arp = y[N_s + Nb:N_s + 2 * Nb]
        CP = y[N_s + 2 * Nb:N_tot]
        Gm, Arpm, CPm = G[:Nm], Arp[:Nm], CP[:Nm]
        LsB = L_s @ B
        D = lambda v: diags(v, 0, format="csc")
        DM = lambda v: coo_matrix((v, (np.arange(Nm), np.arange(Nm))), shape=(Nm, Nb)).tocsc()
        DS = lambda v: coo_matrix((s2b * v, (np.arange(Nm), np.arange(Nm))), shape=(Nb, Nm)).tocsc()
        DBB = lambda v: coo_matrix((s2b * v, (np.arange(Nm), np.arange(Nm))), shape=(Nb, Nb)).tocsc()

        J00 = L_s_V + D(-k1 * Gm)
        J01 = D(np.full(Nm, k2))
        J04 = D(np.full(Nm, kre))
        J05 = DM(-k1 * V)
        J10 = D(k1 * Gm)
        J11 = L_s_Vs + D(-k2 - k3 * Arpm)
        J12 = D(np.full(Nm, k4))
        J15 = DM(k1 * V)
        J16 = DM(-k3 * Vs)
        J21 = D(k3 * Arpm)
        J22 = L_s_A + D(-k4 - kArp * B)
        J23 = D(-kArp * As)
        J26 = DM(k3 * Vs)
        J32 = D(kArp * B)
        J33 = (-DB * L_s.multiply((Gm / G0)[:, None])) + D(kArp * As - kCap * CPm - kDeg)
        J35 = DM(-DB / G0 * LsB)
        J37 = DM(-kCap * B)
        J42 = D(kArp * B)
        J43 = D(kArp * As)
        J44 = D(np.full(Nm, -kre))
        J50 = DS(-k1 * Gm)
        J51 = DS(np.full(Nm, cb * k2))
        J53 = DS(-kon * Gm + koff)
        J55 = L_b_G + DBB(-k1 * V - kon * B)
        J61 = DS(-k3 * Arpm)
        J62 = DS(np.full(Nm, cb * k4))
        J66 = L_b_A + DBB(-k3 * Vs)
        J73 = DS(-kCap * CPm)
        J77 = L_b_C + DBB(-kCap * B)
        return bmat([
            [J00, J01, None, None, J04, J05, None, None],
            [J10, J11, J12, None, None, J15, J16, None],
            [None, J21, J22, J23, None, None, J26, None],
            [None, None, J32, J33, None, J35, None, J37],
            [None, None, J42, J43, J44, None, None, None],
            [J50, J51, None, J53, None, J55, None, None],
            [None, J61, J62, None, None, None, J66, None],
            [None, None, None, J73, None, None, None, J77],
        ], format="csc")

    B0 = B0_norm * np.exp(-r_m**2 / (2.0 * sigma**2))
    y0 = np.concatenate([
        np.full(Nm, V0), np.zeros(Nm), np.zeros(Nm), B0, np.zeros(Nm),
        np.full(Nb, G0), np.full(Nb, Arp0), np.full(Nb, CP0),
    ])
    t_eval = np.linspace(0.0, t_end, max(81, int(round(t_end * 2)) + 1))
    sol = solve_ivp(rhs, [0.0, t_end], y0, method="BDF", jac=jac, rtol=rtol, atol=atol, t_eval=t_eval)
    if not sol.success:
        raise RuntimeError(f"BDF integration failed: {sol.message}")

    B = sol.y[3 * Nm:4 * Nm, :]
    Vtot = sol.y[0:Nm, :] + sol.y[Nm:2 * Nm, :] + sol.y[2 * Nm:3 * Nm, :] + sol.y[4 * Nm:5 * Nm, :]
    B_final = B[:, -1]
    V_final = Vtot[:, -1]
    mask_pk = t_eval >= min(2.0, 0.25 * t_end)
    idx_max_r = int(np.argmax(B_final))
    b_center = float(B_final[0])
    b_max = float(B_final[idx_max_r])
    v_center = float(V_final[0])
    v_edge = float(V_final[-1])
    return {
        "B_peak": float(np.max(B[:, mask_pk])),
        "B_final_center": b_center,
        "B_final_max": b_max,
        "r_peak_final": float(r_m[idx_max_r]),
        "rim_to_center_ratio": float(b_max / b_center) if b_center > 1e-12 else 1.0,
        "V_tot_final_center": v_center,
        "V_tot_final_edge": v_edge,
        "vca_fold_enrichment": float(v_center / v_edge) if v_edge > 1e-12 else 1.0,
    }


def run(inputs, settings):
    return _solve_cone(inputs, settings, conserved_bc=False, surf_axisym=True)

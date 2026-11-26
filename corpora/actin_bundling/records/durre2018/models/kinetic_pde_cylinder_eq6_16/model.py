"""Axisymmetric 2D(r,z) P1-FEM solver for the infinite-pool cylinder kinetic PDEs Eqs. (6)-(16) as printed.

Implements printed boundary conditions Eqs. (14)-(15) verbatim (without +k_2*V* and +k_4*A*).
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import bmat, coo_matrix, diags, eye, kron

NA_PER_UM3_PER_UM = 602.214076


def _solve_cylinder(inputs, settings, conserved_bc, surf_axisym):
    G0 = float(inputs["G_0"])
    Arp0 = float(inputs["Arp_0"])
    CP0 = float(inputs["CP_0"]) * 1e-3  # nM -> uM
    V0 = float(inputs["V_0"])
    B0_norm = float(inputs["B_0_norm"])
    sigma = float(inputs["sigma"]) * 1e-3  # nm -> um
    R = float(inputs["r_cyl"])
    h = float(inputs["h_cyl"])
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

    mesh_n = max(8, int(round(float(settings.get("mesh_n", 24)))))
    Nr = mesh_n
    Nz = max(8, mesh_n // 2)
    r_stretch = float(settings.get("r_stretch", 5.5))
    z_stretch = float(settings.get("z_stretch", 2.5))
    t_end = float(settings.get("t_end", 100.0))
    rtol = float(settings.get("rtol", 1e-5))
    atol = float(settings.get("atol", 1e-7))

    s_r = np.linspace(0.0, 1.0, Nr + 1)
    r = R * np.sinh(r_stretch * s_r) / np.sinh(r_stretch) if r_stretch > 0 else R * s_r
    s_z = np.linspace(0.0, 1.0, Nz + 1)
    z = h * np.sinh(z_stretch * s_z) / np.sinh(z_stretch) if z_stretch > 0 else h * s_z

    Nm = Nr + 1
    dr = np.diff(r)
    rm = 0.5 * (r[:-1] + r[1:])
    i0, i1 = np.arange(Nr), np.arange(1, Nm)

    M_s_ax = np.zeros(Nm)
    np.add.at(M_s_ax, i0, 2.0 * np.pi * dr / 6.0 * (2.0 * r[:-1] + r[1:]))
    np.add.at(M_s_ax, i1, 2.0 * np.pi * dr / 6.0 * (r[:-1] + 2.0 * r[1:]))
    sf_ax = 2.0 * np.pi * rm / dr
    K_s_ax = coo_matrix(
        (np.concatenate([sf_ax, -sf_ax, -sf_ax, sf_ax]),
         (np.concatenate([i0, i0, i1, i1]), np.concatenate([i0, i1, i0, i1]))),
        shape=(Nm, Nm),
    ).tocsc()
    L_b_r = K_s_ax.multiply(1.0 / M_s_ax[:, None]).tocsc()

    if surf_axisym:
        L_s = L_b_r
    else:
        M_s_c = np.zeros(Nm)
        np.add.at(M_s_c, i0, 0.5 * dr)
        np.add.at(M_s_c, i1, 0.5 * dr)
        sf_c = 1.0 / dr
        K_s_c = coo_matrix(
            (np.concatenate([sf_c, -sf_c, -sf_c, sf_c]),
             (np.concatenate([i0, i0, i1, i1]), np.concatenate([i0, i1, i0, i1]))),
            shape=(Nm, Nm),
        ).tocsc()
        L_s = K_s_c.multiply(1.0 / M_s_c[:, None]).tocsc()

    zi = z[:-1]
    Nzi = Nz
    dz = np.diff(zi)
    zf = np.concatenate(([0.0], 0.5 * (zi[:-1] + zi[1:]), [h]))
    dzc = np.diff(zf)
    czp = np.zeros(Nzi)
    czm = np.zeros(Nzi)
    czp[:-1] = (1.0 / dz) / dzc[:-1]
    czm[1:] = (1.0 / dz) / dzc[1:]
    cz_top = (1.0 / (h - zi[-1])) / dzc[-1]
    diag_z = czp + czm
    diag_z[-1] += cz_top
    L_z = diags([-czm[1:], diag_z, -czp[:-1]], [-1, 0, 1], shape=(Nzi, Nzi), format="csc")

    Nb = Nm * Nzi
    L_b = (kron(L_b_r, eye(Nzi, format="csc"), format="csc") + kron(eye(Nm, format="csc"), L_z, format="csc")).tocsc()
    top_idx = np.arange(Nm) * Nzi + (Nzi - 1)
    top_src = np.zeros(Nb)
    top_src[top_idx] = cz_top
    mem_idx = np.arange(Nm) * Nzi
    s2b = (1.0 / dzc[0]) / NA_PER_UM3_PER_UM

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
        Gm, Arpm, CPm = G[mem_idx], Arp[mem_idx], CP[mem_idx]

        dV = L_s_V @ V - k1 * Gm * V + k2 * Vs + kre * Vp
        dVs = L_s_Vs @ Vs + k1 * Gm * V - k2 * Vs - k3 * Arpm * Vs + k4 * As
        dAs = L_s_A @ As + k3 * Arpm * Vs - k4 * As - kArp * As * B
        dB = -DB * (Gm / G0) * (L_s @ B) + kArp * As * B - kCap * CPm * B - kDeg * B
        dVp = kArp * As * B - kre * Vp

        dG = L_b_G @ G + DG * G0 * top_src
        dArp = L_b_A @ Arp + DArp * Arp0 * top_src
        dCP = L_b_C @ CP + DCP * CP0 * top_src
        dG[mem_idx] += s2b * (-k1 * Gm * V + cb * k2 * Vs - kon * Gm * B + koff * B)
        dArp[mem_idx] += s2b * (-k3 * Arpm * Vs + cb * k4 * As)
        dCP[mem_idx] += s2b * (-kCap * CPm * B)
        return np.concatenate([dV, dVs, dAs, dB, dVp, dG, dArp, dCP])

    def jac(t, y):
        V, Vs, As, B, Vp = y[0:Nm], y[Nm:2 * Nm], y[2 * Nm:3 * Nm], y[3 * Nm:4 * Nm], y[4 * Nm:5 * Nm]
        G = y[N_s:N_s + Nb]
        Arp = y[N_s + Nb:N_s + 2 * Nb]
        CP = y[N_s + 2 * Nb:N_tot]
        Gm, Arpm, CPm = G[mem_idx], Arp[mem_idx], CP[mem_idx]
        LsB = L_s @ B
        D = lambda v: diags(v, 0, format="csc")
        DM = lambda v: coo_matrix((v, (np.arange(Nm), mem_idx)), shape=(Nm, Nb)).tocsc()
        DS = lambda v: coo_matrix((s2b * v, (mem_idx, np.arange(Nm))), shape=(Nb, Nm)).tocsc()
        DBB = lambda v: coo_matrix((s2b * v, (mem_idx, mem_idx)), shape=(Nb, Nb)).tocsc()

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

    B0 = B0_norm * np.exp(-r**2 / (2.0 * sigma**2))
    y0 = np.concatenate([
        np.full(Nm, V0), np.zeros(Nm), np.zeros(Nm), B0, np.zeros(Nm),
        np.full(Nb, G0), np.full(Nb, Arp0), np.full(Nb, CP0),
    ])
    t_eval = np.linspace(0.0, t_end, max(201, int(round(t_end * 4)) + 1))
    sol = solve_ivp(rhs, [0.0, t_end], y0, method="BDF", jac=jac, rtol=rtol, atol=atol, t_eval=t_eval)
    if not sol.success:
        raise RuntimeError(f"BDF integration failed: {sol.message}")

    B0_t = sol.y[3 * Nm, :]
    As0_t = sol.y[2 * Nm, :]
    br0_t = kArp * As0_t * B0_t
    mask_fig5 = (t_eval >= 2.0) & (t_eval <= min(40.0, t_end) + 1e-9)
    idx_40 = int(np.argmin(np.abs(t_eval - min(40.0, t_end))))
    return {
        "B_peak": float(np.max(B0_t[mask_fig5])),
        "B_40": float(B0_t[idx_40]),
        "B_100": float(B0_t[-1]),
        "A_star_peak": float(np.max(As0_t[t_eval <= min(40.0, t_end) + 1e-9])),
        "A_star_40": float(As0_t[idx_40]),
        "branching_rate_peak": float(np.max(br0_t[mask_fig5])),
        "branching_rate_40": float(br0_t[idx_40]),
    }


def run(inputs, settings):
    return _solve_cylinder(inputs, settings, conserved_bc=False, surf_axisym=True)

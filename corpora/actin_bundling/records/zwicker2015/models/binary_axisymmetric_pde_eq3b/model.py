import math
import numpy as np
from scipy.fft import dctn, idctn


def run(inputs, settings):
    kf_h = inputs["kf_norm"]
    kb_h = inputs["kb_norm"]
    Vs = inputs["Vs_over_w3"]
    phi_B_bar_in = inputs["phi_B_bar"]
    aspect = inputs["aspect"]
    sep_frac = inputs["sep_frac"]
    dr_frac = inputs["dr_frac"]
    T_tot = inputs["T_k0"]

    nr = int(round(settings.get("nr", 32)))
    nz = int(round(settings.get("nz", 128)))
    dt = float(settings.get("dt", 1.0))

    R_cyl = (Vs / (2.0 * math.pi * aspect)) ** (1.0 / 3.0)
    L_z = 2.0 * aspect * R_cyl
    dr = R_cyl / nr
    dz = L_z / nz
    r = (np.arange(nr) + 0.5) * dr
    z = (np.arange(nz) + 0.5) * dz

    # Conservative axisymmetric radial Laplacian matrix with homogeneous Neumann BCs
    L_r = np.zeros((nr, nr), dtype=float)
    for j in range(nr):
        rp = (j + 1) * dr if j < nr - 1 else 0.0
        rm = j * dr if j > 0 else 0.0
        if j < nr - 1:
            c = rp / (r[j] * dr**2)
            L_r[j, j + 1] = c
            L_r[j, j] -= c
        if j > 0:
            c = rm / (r[j] * dr**2)
            L_r[j, j - 1] = c
            L_r[j, j] -= c

    sqrt_r = np.sqrt(r)
    S_r = (sqrt_r[:, None] * L_r) / sqrt_r[None, :]
    ev_r, Q_r = np.linalg.eigh(S_r)
    ev_z = -2.0 * (1.0 - np.cos(np.pi * np.arange(nz) / nz)) / (dz**2)
    lam_lap = ev_r[:, None] + ev_z[None, :]

    def fwd(u):
        return dctn(Q_r.T @ (sqrt_r[:, None] * u), type=2, axes=(1,), norm="ortho")

    def inv(uh):
        return (Q_r @ idctn(uh, type=2, axes=(1,), norm="ortho")) / sqrt_r[:, None]

    passive = (kf_h == 0.0 and kb_h == 0.0)
    phi_bar = phi_B_bar_in if passive else (kf_h / (kf_h + kb_h))
    Vbar_est = Vs * phi_bar / 2.0
    Rbar_est = (3.0 * Vbar_est / (4.0 * math.pi)) ** (1.0 / 3.0)
    R1 = Rbar_est * (1.0 + dr_frac)
    R2 = Rbar_est * (1.0 - dr_frac)
    z1 = 0.5 * L_z * (1.0 - sep_frac)
    z2 = 0.5 * L_z * (1.0 + sep_frac)

    RR, ZZ = np.meshgrid(r, z, indexing="ij")
    d1 = np.sqrt(RR**2 + (ZZ - z1) ** 2)
    d2 = np.sqrt(RR**2 + (ZZ - z2) ** 2)
    phi = 0.5 * (1.0 - np.tanh((d1 - R1) / 1.0)) + 0.5 * (1.0 - np.tanh((d2 - R2) / 1.0))
    phi = np.clip(phi, 0.0, 1.0)
    vol_w = 2.0 * math.pi * r[:, None] * dr * dz

    if passive:
        curr = float(np.sum(phi * vol_w) / Vs)
        if curr < phi_bar:
            phi += (phi_bar - curr) * (1.0 - phi) / float(np.sum((1.0 - phi) * vol_w) / Vs)

    # Second-order IMEX (Crank-Nicolson linear + Heun/trapezoidal nonlinear)
    S_stab = 2.0
    L_op = S_stab * lam_lap - 0.25 * (lam_lap**2) - (kf_h + kb_h)
    nsteps = max(1, int(round(T_tot / dt)))
    h = T_tot / nsteps
    lhs_e = 1.0 - h * L_op
    lhs_cn = 1.0 - 0.5 * h * L_op
    rhs_cn = 1.0 + 0.5 * h * L_op
    kf_hat = fwd(np.full((nr, nz), kf_h, dtype=float))

    def nl_hat_fn(p):
        nl = (1.0 - S_stab) * p - 3.0 * (p**2) + 2.0 * (p**3)
        return lam_lap * fwd(nl) + kf_hat

    phi_hat = fwd(phi)
    sample_every = max(1, nsteps // 30)
    ts, dVs, Vmeans = [], [], []
    phi_min = float(np.min(phi))

    for step in range(nsteps + 1):
        if step % sample_every == 0:
            ind = (phi**2) * (3.0 - 2.0 * phi)
            V1 = float(np.sum(ind[:, : nz // 2] * vol_w))
            V2 = float(np.sum(ind[:, nz // 2 :] * vol_w))
            ts.append(step * h)
            dVs.append(abs(V1 - V2))
            Vmeans.append(0.5 * (V1 + V2))
            pmin = float(np.min(phi))
            if pmin < phi_min:
                phi_min = pmin
        if step < nsteps:
            nh0 = nl_hat_fn(phi)
            phi_pred = inv((phi_hat + h * nh0) / lhs_e)
            nh1 = nl_hat_fn(phi_pred)
            phi_hat = (rhs_cn * phi_hat + 0.5 * h * (nh0 + nh1)) / lhs_cn
            phi = inv(phi_hat)

    ts_arr = np.array(ts)
    dVs_arr = np.array(dVs)
    Vm_arr = np.array(Vmeans)
    t_fit_start = T_tot / 3.0
    mask = (ts_arr >= t_fit_start) & (dVs_arr > 1e-6)
    if np.sum(mask) < 3:
        return {"lam_norm": None, "V_bar_over_w3": float(np.mean(Vm_arr)), "phi_min": phi_min}

    slope = float(np.polyfit(ts_arr[mask], np.log(dVs_arr[mask]), 1)[0])
    V_bar = float(np.mean(Vm_arr[mask]))
    return {
        "lam_norm": slope,
        "V_bar_over_w3": V_bar,
        "phi_min": phi_min,
    }

"""2nd-order Strang-splitting variant (Brownian dt/2 -> deterministic spring dt -> Brownian dt/2)
corresponding to maxian2022 Sec. 2.4 (p. 5), contrasted with TemporalIntegrator.py:441-499's
1st-order Lie-Trotter implementation. Outputs both midpoint (after step 1, per Fig. S1 caption)
and true full-step endpoint (after step 3) spring extension statistics."""

import math
import numpy as np


def _cheb_kind2(Nx, L):
    th = np.flipud(np.pi * np.arange(Nx)) / (Nx - 1)
    c = 2.0 / np.concatenate(([1.0], 1.0 - np.arange(2, Nx, 2) ** 2))
    c = np.concatenate((c, c[Nx // 2 - 1 : 0 : -1]))
    w = np.real(np.fft.ifft(c))
    w[0] /= 2.0
    w = np.concatenate((w, [w[0]]))
    s = 0.5 * L * (np.cos(th) + 1.0)
    w = 0.5 * L * w
    return s, w


def _regularize_s(s0, delta, L):
    sNew = 0.5 * L * np.ones_like(s0)
    if 0.0 < delta < 0.5:
        x = 2.0 * s0 / L - 1.0
        regwt = np.tanh((x + 1.0) / delta) - np.tanh((x - 1.0) / delta) - 1.0
        sNew = s0.copy()
        mask = s0 < 0.5 * L
        sNew[mask] = regwt[mask] * s0[mask] + (1.0 - regwt[mask] ** 2) * delta * L / 2.0
        sNew[s0 > 0.5 * L] = L - np.flip(sNew[mask])
    return sNew


def _local_drag_mobilities(L, eps, mu, delta, sigma_rel, N, smooth_force):
    Nx = N + 1
    s, w = _cheb_kind2(Nx, L)
    sbar = _regularize_s(s, delta, L)
    a = eps * L
    c_val = np.log(4.0 * sbar * (L - sbar) / (a * a))
    m_perp = (c_val + 1.0) / (8.0 * np.pi * mu)
    m_par = (2.0 * c_val - 2.0) / (8.0 * np.pi * mu)
    sc = s - 0.5 * L
    N_perp = 1.0 / float(np.dot(w, 1.0 / m_perp))
    N_par = 1.0 / float(np.dot(w, 1.0 / m_par))
    N_rot = 1.0 / float(np.dot(w, (sc ** 2) / m_perp))
    if smooth_force:
        sig = sigma_rel * L
        dh = np.exp(-((s - L) ** 2) / (2.0 * sig * sig))
        dh_norm = dh / float(np.dot(w, dh))
        arm = float(np.dot(w, sc * dh_norm))
    else:
        arm = 0.5 * L
    return N_perp, N_par, N_rot, arm


def run(inputs, settings):
    L = float(inputs["L"])
    eps = float(inputs["eps"])
    mu = float(inputs["mu"])
    K_c = float(inputs["K_c"])
    ell = float(inputs["ell"])
    kBT = float(inputs["kBT"])
    delta = float(inputs["delta"])
    sigma_rel = float(inputs["sigma_rel"])

    dt = float(settings.get("dt", 0.0005))
    N = int(round(float(settings.get("N", 50))))
    smooth_force = bool(settings.get("smooth_force", True))
    t_sim = float(settings.get("t_sim", 40.0))
    burn_t = float(settings.get("burn_t", 1.0))
    seed = int(settings.get("seed", 0))

    N_perp, N_par, N_rot, arm = _local_drag_mobilities(
        L, eps, mu, delta, sigma_rel, N, smooth_force
    )
    dN_t = N_par - N_perp

    # Half-step Brownian amplitudes (dt / 2)
    half_dt = 0.5 * dt
    amp_perp = math.sqrt(2.0 * kBT * N_perp * half_dt)
    amp_par = math.sqrt(2.0 * kBT * N_par * half_dt)
    d_amp_t = amp_par - amp_perp
    amp_rot = math.sqrt(2.0 * kBT * N_rot * half_dt)

    n_steps = max(1, int(round(t_sim / dt)))
    burn_steps = max(0, int(round(burn_t / dt)))

    rng = np.random.default_rng(seed)
    W_all = rng.standard_normal((n_steps, 24)).tolist()

    m1x = 0.0; m1y = 0.0; m1z = 0.0
    m2x = 0.0; m2y = ell; m2z = 0.0
    t1x = 1.0; t1y = 0.0; t1z = 0.0
    t2x = 1.0; t2y = 0.0; t2z = 0.0

    s_end = 0.0; s2_end = 0.0
    s_bd = 0.0; s2_bd = 0.0
    n_samp = 0

    sqrt = math.sqrt
    cos = math.cos
    sin = math.sin
    n_rot_dt = N_rot * dt

    for step in range(n_steps):
        w = W_all[step]

        # --- Step 1: BrownianUpdate(dt/2) ---
        dt1 = w[3] * t1x + w[4] * t1y + w[5] * t1z
        m1x += amp_perp * w[3] + d_amp_t * dt1 * t1x
        m1y += amp_perp * w[4] + d_amp_t * dt1 * t1y
        m1z += amp_perp * w[5] + d_amp_t * dt1 * t1z
        dot1 = w[0] * t1x + w[1] * t1y + w[2] * t1z
        ox = amp_rot * (w[0] - dot1 * t1x)
        oy = amp_rot * (w[1] - dot1 * t1y)
        oz = amp_rot * (w[2] - dot1 * t1z)
        th = sqrt(ox * ox + oy * oy + oz * oz)
        if th > 1e-12:
            inv_th = 1.0 / th
            ux = ox * inv_th; uy = oy * inv_th; uz = oz * inv_th
            ct = cos(th); st = sin(th)
            cx = uy * t1z - uz * t1y
            cy = uz * t1x - ux * t1z
            cz = ux * t1y - uy * t1x
            t1x = t1x * ct + cx * st
            t1y = t1y * ct + cy * st
            t1z = t1z * ct + cz * st
            inv_n = 1.0 / sqrt(t1x * t1x + t1y * t1y + t1z * t1z)
            t1x *= inv_n; t1y *= inv_n; t1z *= inv_n

        dt2 = w[9] * t2x + w[10] * t2y + w[11] * t2z
        m2x += amp_perp * w[9] + d_amp_t * dt2 * t2x
        m2y += amp_perp * w[10] + d_amp_t * dt2 * t2y
        m2z += amp_perp * w[11] + d_amp_t * dt2 * t2z
        dot2 = w[6] * t2x + w[7] * t2y + w[8] * t2z
        ox = amp_rot * (w[6] - dot2 * t2x)
        oy = amp_rot * (w[7] - dot2 * t2y)
        oz = amp_rot * (w[8] - dot2 * t2z)
        th = sqrt(ox * ox + oy * oy + oz * oz)
        if th > 1e-12:
            inv_th = 1.0 / th
            ux = ox * inv_th; uy = oy * inv_th; uz = oz * inv_th
            ct = cos(th); st = sin(th)
            cx = uy * t2z - uz * t2y
            cy = uz * t2x - ux * t2z
            cz = ux * t2y - uy * t2x
            t2x = t2x * ct + cx * st
            t2y = t2y * ct + cy * st
            t2z = t2z * ct + cz * st
            inv_n = 1.0 / sqrt(t2x * t2x + t2y * t2y + t2z * t2z)
            t2x *= inv_n; t2y *= inv_n; t2z *= inv_n

        # Midpoint separation (after step 1, Fig. S1 caption)
        dx = (m2x + arm * t2x) - (m1x + arm * t1x)
        dy = (m2y + arm * t2y) - (m1y + arm * t1y)
        dz = (m2z + arm * t2z) - (m1z + arm * t1z)
        r_mid = sqrt(dx * dx + dy * dy + dz * dz)

        # --- Step 2: Deterministic explicit spring update(dt) ---
        if r_mid > 1e-15:
            fac = K_c * (1.0 - ell / r_mid)
            fx = fac * dx; fy = fac * dy; fz = fac * dz
        else:
            fx = 0.0; fy = 0.0; fz = 0.0

        fdot1 = fx * t1x + fy * t1y + fz * t1z
        m1x += dt * (N_perp * fx + dN_t * fdot1 * t1x)
        m1y += dt * (N_perp * fy + dN_t * fdot1 * t1y)
        m1z += dt * (N_perp * fz + dN_t * fdot1 * t1z)
        ox = n_rot_dt * arm * (t1y * fz - t1z * fy)
        oy = n_rot_dt * arm * (t1z * fx - t1x * fz)
        oz = n_rot_dt * arm * (t1x * fy - t1y * fx)
        th = sqrt(ox * ox + oy * oy + oz * oz)
        if th > 1e-12:
            inv_th = 1.0 / th
            ux = ox * inv_th; uy = oy * inv_th; uz = oz * inv_th
            ct = cos(th); st = sin(th)
            cx = uy * t1z - uz * t1y
            cy = uz * t1x - ux * t1z
            cz = ux * t1y - uy * t1x
            t1x = t1x * ct + cx * st
            t1y = t1y * ct + cy * st
            t1z = t1z * ct + cz * st
            inv_n = 1.0 / sqrt(t1x * t1x + t1y * t1y + t1z * t1z)
            t1x *= inv_n; t1y *= inv_n; t1z *= inv_n

        fdot2 = fx * t2x + fy * t2y + fz * t2z
        m2x -= dt * (N_perp * fx + dN_t * fdot2 * t2x)
        m2y -= dt * (N_perp * fy + dN_t * fdot2 * t2y)
        m2z -= dt * (N_perp * fz + dN_t * fdot2 * t2z)
        ox = -n_rot_dt * arm * (t2y * fz - t2z * fy)
        oy = -n_rot_dt * arm * (t2z * fx - t2x * fz)
        oz = -n_rot_dt * arm * (t2x * fy - t2y * fx)
        th = sqrt(ox * ox + oy * oy + oz * oz)
        if th > 1e-12:
            inv_th = 1.0 / th
            ux = ox * inv_th; uy = oy * inv_th; uz = oz * inv_th
            ct = cos(th); st = sin(th)
            cx = uy * t2z - uz * t2y
            cy = uz * t2x - ux * t2z
            cz = ux * t2y - uy * t2x
            t2x = t2x * ct + cx * st
            t2y = t2y * ct + cy * st
            t2z = t2z * ct + cz * st
            inv_n = 1.0 / sqrt(t2x * t2x + t2y * t2y + t2z * t2z)
            t2x *= inv_n; t2y *= inv_n; t2z *= inv_n

        # --- Step 3: BrownianUpdate(dt/2) ---
        dt1 = w[15] * t1x + w[16] * t1y + w[17] * t1z
        m1x += amp_perp * w[15] + d_amp_t * dt1 * t1x
        m1y += amp_perp * w[16] + d_amp_t * dt1 * t1y
        m1z += amp_perp * w[17] + d_amp_t * dt1 * t1z
        dot1 = w[12] * t1x + w[13] * t1y + w[14] * t1z
        ox = amp_rot * (w[12] - dot1 * t1x)
        oy = amp_rot * (w[13] - dot1 * t1y)
        oz = amp_rot * (w[14] - dot1 * t1z)
        th = sqrt(ox * ox + oy * oy + oz * oz)
        if th > 1e-12:
            inv_th = 1.0 / th
            ux = ox * inv_th; uy = oy * inv_th; uz = oz * inv_th
            ct = cos(th); st = sin(th)
            cx = uy * t1z - uz * t1y
            cy = uz * t1x - ux * t1z
            cz = ux * t1y - uy * t1x
            t1x = t1x * ct + cx * st
            t1y = t1y * ct + cy * st
            t1z = t1z * ct + cz * st
            inv_n = 1.0 / sqrt(t1x * t1x + t1y * t1y + t1z * t1z)
            t1x *= inv_n; t1y *= inv_n; t1z *= inv_n

        dt2 = w[21] * t2x + w[22] * t2y + w[23] * t2z
        m2x += amp_perp * w[21] + d_amp_t * dt2 * t2x
        m2y += amp_perp * w[22] + d_amp_t * dt2 * t2y
        m2z += amp_perp * w[23] + d_amp_t * dt2 * t2z
        dot2 = w[18] * t2x + w[19] * t2y + w[20] * t2z
        ox = amp_rot * (w[18] - dot2 * t2x)
        oy = amp_rot * (w[19] - dot2 * t2y)
        oz = amp_rot * (w[20] - dot2 * t2z)
        th = sqrt(ox * ox + oy * oy + oz * oz)
        if th > 1e-12:
            inv_th = 1.0 / th
            ux = ox * inv_th; uy = oy * inv_th; uz = oz * inv_th
            ct = cos(th); st = sin(th)
            cx = uy * t2z - uz * t2y
            cy = uz * t2x - ux * t2z
            cz = ux * t2y - uy * t2x
            t2x = t2x * ct + cx * st
            t2y = t2y * ct + cy * st
            t2z = t2z * ct + cz * st
            inv_n = 1.0 / sqrt(t2x * t2x + t2y * t2y + t2z * t2z)
            t2x *= inv_n; t2y *= inv_n; t2z *= inv_n

        if step >= burn_steps:
            dx = (m2x + arm * t2x) - (m1x + arm * t1x)
            dy = (m2y + arm * t2y) - (m1y + arm * t1y)
            dz = (m2z + arm * t2z) - (m1z + arm * t1z)
            r_end = sqrt(dx * dx + dy * dy + dz * dz)
            s_bd += r_mid
            s2_bd += r_mid * r_mid
            s_end += r_end
            s2_end += r_end * r_end
            n_samp += 1

    mean_bd = s_bd / n_samp
    var_bd = max(0.0, s2_bd / n_samp - mean_bd * mean_bd)
    mean_end = s_end / n_samp
    var_end = max(0.0, s2_end / n_samp - mean_end * mean_end)

    return {
        "mean_x_end": float(mean_end),
        "std_x_end": float(math.sqrt(var_end)),
        "mean_x_after_bd": float(mean_bd),
        "std_x_after_bd": float(math.sqrt(var_bd)),
    }

import math
import random


def run(inputs, settings):
    alpha = float(inputs["alpha"])
    beta = float(inputs["beta"])
    gamma = float(inputs["gamma"])
    Kc = float(inputs["Kc"])
    ell = float(inputs["ell"])
    mu = float(inputs["mu"])
    L = float(inputs["L"])
    kBT = float(inputs["kBT"])

    dt = float(settings.get("dt", 0.0025))
    t_max = float(settings.get("t_max", 100.0))
    burn_t = float(settings.get("burn_t", 1.0))
    seed = int(settings.get("seed", 0))

    rng = random.Random(seed)
    gauss = rng.gauss
    sqrt = math.sqrt
    sin = math.sin
    cos = math.cos

    # Initial state (SM Sec. C): X^(1)(s) = (s, 0, 0), X^(2)(s) = (s, ell, 0)
    half_L = 0.5 * L
    x0, y0, z0 = half_L, 0.0, 0.0
    x1, y1, z1 = half_L, ell, 0.0
    tx0, ty0, tz0 = 1.0, 0.0, 0.0
    tx1, ty1, tz1 = 1.0, 0.0, 0.0

    # 1 Pa*s = 1 pN*s/um^2
    inv_muL = 1.0 / (mu * L)
    inv_muL3 = 1.0 / (mu * L * L * L)
    s_a = sqrt(alpha)
    s_diff = sqrt(alpha + beta) - s_a
    s_g = sqrt(gamma)
    half_dt = 0.5 * dt
    b_t = sqrt(2.0 * kBT * half_dt * inv_muL)
    b_r = sqrt(2.0 * kBT * half_dt * inv_muL3) * s_g

    n_steps = int(round(t_max / dt))
    burn = int(round(burn_t / dt))
    sum_x = 0.0
    sum_x2 = 0.0
    count = 0

    for step in range(n_steps):
        # Step 1: Brownian dt/2 (Eq. 15 & Eq. S12) for fiber 0
        w0, w1, w2 = gauss(0.0, 1.0), gauss(0.0, 1.0), gauss(0.0, 1.0)
        dot = tx0 * w0 + ty0 * w1 + tz0 * w2
        x0 += b_t * (s_a * w0 + s_diff * dot * tx0)
        y0 += b_t * (s_a * w1 + s_diff * dot * ty0)
        z0 += b_t * (s_a * w2 + s_diff * dot * tz0)
        w0, w1, w2 = gauss(0.0, 1.0), gauss(0.0, 1.0), gauss(0.0, 1.0)
        dot = tx0 * w0 + ty0 * w1 + tz0 * w2
        rx, ry, rz = b_r * (w0 - dot * tx0), b_r * (w1 - dot * ty0), b_r * (w2 - dot * tz0)
        ang = sqrt(rx * rx + ry * ry + rz * rz)
        if ang > 1e-15:
            ia = 1.0 / ang
            kx, ky, kz = rx * ia, ry * ia, rz * ia
            ca, sa = cos(ang), sin(ang)
            omc = 1.0 - ca
            kdt = kx * tx0 + ky * ty0 + kz * tz0
            nx = tx0 * ca + (ky * tz0 - kz * ty0) * sa + kx * kdt * omc
            ny = ty0 * ca + (kz * tx0 - kx * tz0) * sa + ky * kdt * omc
            nz = tz0 * ca + (kx * ty0 - ky * tx0) * sa + kz * kdt * omc
            inorm = 1.0 / sqrt(nx * nx + ny * ny + nz * nz)
            tx0, ty0, tz0 = nx * inorm, ny * inorm, nz * inorm

        # Step 1: Brownian dt/2 for fiber 1
        w0, w1, w2 = gauss(0.0, 1.0), gauss(0.0, 1.0), gauss(0.0, 1.0)
        dot = tx1 * w0 + ty1 * w1 + tz1 * w2
        x1 += b_t * (s_a * w0 + s_diff * dot * tx1)
        y1 += b_t * (s_a * w1 + s_diff * dot * ty1)
        z1 += b_t * (s_a * w2 + s_diff * dot * tz1)
        w0, w1, w2 = gauss(0.0, 1.0), gauss(0.0, 1.0), gauss(0.0, 1.0)
        dot = tx1 * w0 + ty1 * w1 + tz1 * w2
        rx, ry, rz = b_r * (w0 - dot * tx1), b_r * (w1 - dot * ty1), b_r * (w2 - dot * tz1)
        ang = sqrt(rx * rx + ry * ry + rz * rz)
        if ang > 1e-15:
            ia = 1.0 / ang
            kx, ky, kz = rx * ia, ry * ia, rz * ia
            ca, sa = cos(ang), sin(ang)
            omc = 1.0 - ca
            kdt = kx * tx1 + ky * ty1 + kz * tz1
            nx = tx1 * ca + (ky * tz1 - kz * ty1) * sa + kx * kdt * omc
            ny = ty1 * ca + (kz * tx1 - kx * tz1) * sa + ky * kdt * omc
            nz = tz1 * ca + (kx * ty1 - ky * tx1) * sa + kz * kdt * omc
            inorm = 1.0 / sqrt(nx * nx + ny * ny + nz * nz)
            tx1, ty1, tz1 = nx * inorm, ny * inorm, nz * inorm

        # Midpoint spring extension measurement (Fig. S1 caption)
        dx = (x0 + half_L * tx0) - (x1 + half_L * tx1)
        dy = (y0 + half_L * ty0) - (y1 + half_L * ty1)
        dz = (z0 + half_L * tz0) - (z1 + half_L * tz1)
        dist = sqrt(dx * dx + dy * dy + dz * dz)
        if step >= burn:
            sum_x += dist
            sum_x2 += dist * dist
            count += 1

        # Step 2: Explicit deterministic spring update over dt
        if dist > 1e-15:
            fac = -Kc * (1.0 - ell / dist)
            fx, fy, fz = fac * dx, fac * dy, fac * dz
        else:
            fx = fy = fz = 0.0

        dot = tx0 * fx + ty0 * fy + tz0 * fz
        x0 += dt * inv_muL * (alpha * fx + beta * dot * tx0)
        y0 += dt * inv_muL * (alpha * fy + beta * dot * ty0)
        z0 += dt * inv_muL * (alpha * fz + beta * dot * tz0)
        qx = half_L * (ty0 * fz - tz0 * fy)
        qy = half_L * (tz0 * fx - tx0 * fz)
        qz = half_L * (tx0 * fy - ty0 * fx)
        rx, ry, rz = dt * inv_muL3 * gamma * qx, dt * inv_muL3 * gamma * qy, dt * inv_muL3 * gamma * qz
        ang = sqrt(rx * rx + ry * ry + rz * rz)
        if ang > 1e-15:
            ia = 1.0 / ang
            kx, ky, kz = rx * ia, ry * ia, rz * ia
            ca, sa = cos(ang), sin(ang)
            omc = 1.0 - ca
            kdt = kx * tx0 + ky * ty0 + kz * tz0
            nx = tx0 * ca + (ky * tz0 - kz * ty0) * sa + kx * kdt * omc
            ny = ty0 * ca + (kz * tx0 - kx * tz0) * sa + ky * kdt * omc
            nz = tz0 * ca + (kx * ty0 - ky * tx0) * sa + kz * kdt * omc
            inorm = 1.0 / sqrt(nx * nx + ny * ny + nz * nz)
            tx0, ty0, tz0 = nx * inorm, ny * inorm, nz * inorm

        fx, fy, fz = -fx, -fy, -fz
        dot = tx1 * fx + ty1 * fy + tz1 * fz
        x1 += dt * inv_muL * (alpha * fx + beta * dot * tx1)
        y1 += dt * inv_muL * (alpha * fy + beta * dot * ty1)
        z1 += dt * inv_muL * (alpha * fz + beta * dot * tz1)
        qx = half_L * (ty1 * fz - tz1 * fy)
        qy = half_L * (tz1 * fx - tx1 * fz)
        qz = half_L * (tx1 * fy - ty1 * fx)
        rx, ry, rz = dt * inv_muL3 * gamma * qx, dt * inv_muL3 * gamma * qy, dt * inv_muL3 * gamma * qz
        ang = sqrt(rx * rx + ry * ry + rz * rz)
        if ang > 1e-15:
            ia = 1.0 / ang
            kx, ky, kz = rx * ia, ry * ia, rz * ia
            ca, sa = cos(ang), sin(ang)
            omc = 1.0 - ca
            kdt = kx * tx1 + ky * ty1 + kz * tz1
            nx = tx1 * ca + (ky * tz1 - kz * ty1) * sa + kx * kdt * omc
            ny = ty1 * ca + (kz * tx1 - kx * tz1) * sa + ky * kdt * omc
            nz = tz1 * ca + (kx * ty1 - ky * tx1) * sa + kz * kdt * omc
            inorm = 1.0 / sqrt(nx * nx + ny * ny + nz * nz)
            tx1, ty1, tz1 = nx * inorm, ny * inorm, nz * inorm

        # Step 3: Brownian dt/2 for fiber 0
        w0, w1, w2 = gauss(0.0, 1.0), gauss(0.0, 1.0), gauss(0.0, 1.0)
        dot = tx0 * w0 + ty0 * w1 + tz0 * w2
        x0 += b_t * (s_a * w0 + s_diff * dot * tx0)
        y0 += b_t * (s_a * w1 + s_diff * dot * ty0)
        z0 += b_t * (s_a * w2 + s_diff * dot * tz0)
        w0, w1, w2 = gauss(0.0, 1.0), gauss(0.0, 1.0), gauss(0.0, 1.0)
        dot = tx0 * w0 + ty0 * w1 + tz0 * w2
        rx, ry, rz = b_r * (w0 - dot * tx0), b_r * (w1 - dot * ty0), b_r * (w2 - dot * tz0)
        ang = sqrt(rx * rx + ry * ry + rz * rz)
        if ang > 1e-15:
            ia = 1.0 / ang
            kx, ky, kz = rx * ia, ry * ia, rz * ia
            ca, sa = cos(ang), sin(ang)
            omc = 1.0 - ca
            kdt = kx * tx0 + ky * ty0 + kz * tz0
            nx = tx0 * ca + (ky * tz0 - kz * ty0) * sa + kx * kdt * omc
            ny = ty0 * ca + (kz * tx0 - kx * tz0) * sa + ky * kdt * omc
            nz = tz0 * ca + (kx * ty0 - ky * tx0) * sa + kz * kdt * omc
            inorm = 1.0 / sqrt(nx * nx + ny * ny + nz * nz)
            tx0, ty0, tz0 = nx * inorm, ny * inorm, nz * inorm

        # Step 3: Brownian dt/2 for fiber 1
        w0, w1, w2 = gauss(0.0, 1.0), gauss(0.0, 1.0), gauss(0.0, 1.0)
        dot = tx1 * w0 + ty1 * w1 + tz1 * w2
        x1 += b_t * (s_a * w0 + s_diff * dot * tx1)
        y1 += b_t * (s_a * w1 + s_diff * dot * ty1)
        z1 += b_t * (s_a * w2 + s_diff * dot * tz1)
        w0, w1, w2 = gauss(0.0, 1.0), gauss(0.0, 1.0), gauss(0.0, 1.0)
        dot = tx1 * w0 + ty1 * w1 + tz1 * w2
        rx, ry, rz = b_r * (w0 - dot * tx1), b_r * (w1 - dot * ty1), b_r * (w2 - dot * tz1)
        ang = sqrt(rx * rx + ry * ry + rz * rz)
        if ang > 1e-15:
            ia = 1.0 / ang
            kx, ky, kz = rx * ia, ry * ia, rz * ia
            ca, sa = cos(ang), sin(ang)
            omc = 1.0 - ca
            kdt = kx * tx1 + ky * ty1 + kz * tz1
            nx = tx1 * ca + (ky * tz1 - kz * ty1) * sa + kx * kdt * omc
            ny = ty1 * ca + (kz * tx1 - kx * tz1) * sa + ky * kdt * omc
            nz = tz1 * ca + (kx * ty1 - ky * tx1) * sa + kz * kdt * omc
            inorm = 1.0 / sqrt(nx * nx + ny * ny + nz * nz)
            tx1, ty1, tz1 = nx * inorm, ny * inorm, nz * inorm

    mean_x = sum_x / count
    std_x = sqrt(max(0.0, sum_x2 / count - mean_x * mean_x))
    return {"mean_x": mean_x, "std_x": std_x}

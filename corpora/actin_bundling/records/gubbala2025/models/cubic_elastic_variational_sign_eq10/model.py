import numpy as np


def run(inputs, settings):
    E3 = float(inputs["E3"])
    N = int(settings.get("n_grid", 32))
    L = float(settings.get("L_box", 32.0))
    eps = float(settings.get("eps_fd", 1e-6))

    dx = L / N
    x = np.arange(N) * dx
    X, Y = np.meshgrid(x, x, indexing="ij")
    kx = 2.0 * np.pi * np.fft.fftfreq(N, d=dx)
    ky = 2.0 * np.pi * np.fft.rfftfreq(N, d=dx)
    KX, KY = np.meshgrid(kx, ky, indexing="ij")

    # Smooth periodic 2D test field
    Qxx = 0.4 * np.sin(2.0 * np.pi * X / L) * np.cos(2.0 * np.pi * Y / L) + 0.3 * np.cos(4.0 * np.pi * X / L)
    Qxy = 0.3 * np.cos(2.0 * np.pi * X / L) * np.sin(2.0 * np.pi * Y / L) - 0.2 * np.sin(2.0 * np.pi * Y / L)

    def eval_F_E3(qxx, qxy):
        qxx_k = np.fft.rfft2(qxx)
        qxy_k = np.fft.rfft2(qxy)
        dxx_x = np.fft.irfft2(1j * KX * qxx_k, s=(N, N))
        dxx_y = np.fft.irfft2(1j * KY * qxx_k, s=(N, N))
        dxy_x = np.fft.irfft2(1j * KX * qxy_k, s=(N, N))
        dxy_y = np.fft.irfft2(1j * KY * qxy_k, s=(N, N))
        gxx = 2.0 * (dxx_x**2 + dxy_x**2)
        gyy = 2.0 * (dxx_y**2 + dxy_y**2)
        gxy = 2.0 * (dxx_x * dxx_y + dxy_x * dxy_y)
        return float(np.sum(0.5 * E3 * (qxx * (gxx - gyy) + 2.0 * qxy * gxy)) * dx * dx)

    def rhs_E3(qxx, qxy, e3_sign):
        qxx_k = np.fft.rfft2(qxx)
        qxy_k = np.fft.rfft2(qxy)
        dxx_x = np.fft.irfft2(1j * KX * qxx_k, s=(N, N))
        dxx_y = np.fft.irfft2(1j * KY * qxx_k, s=(N, N))
        dxy_x = np.fft.irfft2(1j * KX * qxy_k, s=(N, N))
        dxy_y = np.fft.irfft2(1j * KY * qxy_k, s=(N, N))
        d2xx_xx = np.fft.irfft2(-KX**2 * qxx_k, s=(N, N))
        d2xx_yy = np.fft.irfft2(-KY**2 * qxx_k, s=(N, N))
        d2xx_xy = np.fft.irfft2(-KX * KY * qxx_k, s=(N, N))
        d2xy_xx = np.fft.irfft2(-KX**2 * qxy_k, s=(N, N))
        d2xy_yy = np.fft.irfft2(-KY**2 * qxy_k, s=(N, N))
        d2xy_xy = np.fft.irfft2(-KX * KY * qxy_k, s=(N, N))
        # Printed Eq. (10) E_3 term corresponds to e3_sign = +1; variational -delta F_E3/delta Q is e3_sign = -1
        rxx = -e3_sign * E3 * (
            qxx * (d2xx_xx - d2xx_yy)
            + 2.0 * qxy * d2xx_xy
            + 0.5 * dxx_x**2
            + dxx_x * dxy_y
            - 0.5 * dxx_y**2
            + dxx_y * dxy_x
            - 0.5 * dxy_x**2
            + 0.5 * dxy_y**2
        )
        rxy = -e3_sign * E3 * (
            qxx * (d2xy_xx - d2xy_yy)
            + 2.0 * qxy * d2xy_xy
            - dxx_x * dxx_y
            + dxx_x * dxy_x
            - dxx_y * dxy_y
            + dxy_x * dxy_y
        )
        return rxx, rxy

    F0 = eval_F_E3(Qxx, Qxy)
    out = {}
    for key, sgn in (("dissipation_ratio_printed", +1.0), ("dissipation_ratio_variational", -1.0)):
        rxx, rxy = rhs_E3(Qxx, Qxy, sgn)
        dF_dt = (eval_F_E3(Qxx + eps * rxx, Qxy + eps * rxy) - F0) / eps
        expected = -2.0 * float(np.sum(rxx**2 + rxy**2) * dx * dx)
        out[key] = float(dF_dt / expected) if expected != 0.0 else float("nan")
    return out

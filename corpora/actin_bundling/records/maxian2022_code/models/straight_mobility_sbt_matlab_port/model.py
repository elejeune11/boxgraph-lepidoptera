"""Faithful Python port of Matlab/MainPrograms/StraightMobility.m (lines 1-45),
RegularizeS.m (lines 5-20), AsymptoticLocalMobility.m (lines 8-62, fixing the broken
getGrandMloc call on StraightMobility.m:21), precomputeStokesletInts.m (lines 7-50),
and StokesletFinitePartMatrix.m (lines 6-39)."""

import numpy as np


def _cheb_pts_wts(N, dom, kind=1):
    a_d, b_d = float(dom[0]), float(dom[1])
    if kind == 1:
        th = np.flipud(2.0 * np.arange(N) + 1.0) * np.pi / (2.0 * N)
        m = 2.0 / np.concatenate(([1.0], 1.0 - np.arange(2, N, 2) ** 2))
        if N % 2 == 1:
            c = np.concatenate((m, -m[(N + 1) // 2 - 1 : 0 : -1]))
        else:
            c = np.concatenate((m, [0.0], -m[N // 2 - 1 : 0 : -1]))
        v = np.exp(1j * np.arange(N) * np.pi / N)
        w = np.real(np.fft.ifft(c * v))
    elif kind == 2:
        th = np.flipud(np.pi * np.arange(N)) / (N - 1)
        c = 2.0 / np.concatenate(([1.0], 1.0 - np.arange(2, N, 2) ** 2))
        c = np.concatenate((c, c[N // 2 - 1 : 0 : -1]))
        w = np.real(np.fft.ifft(c))
        w[0] /= 2.0
        w = np.concatenate((w, [w[0]]))
    else:
        raise ValueError(f"Unsupported Chebyshev kind: {kind}")
    x = 0.5 * (b_d - a_d) * (np.cos(th) + 1.0) + a_d
    w = 0.5 * (b_d - a_d) * w
    return x, w, th


def _diff_mat(N, dom, kind=1):
    _, _, th = _cheb_pts_wts(N, dom, kind)
    CtoV = np.cos(np.outer(th, np.arange(N)))
    VtoC = np.linalg.inv(CtoV)
    Dcos = np.zeros_like(VtoC)
    Dcos[N - 2, :] = 2.0 * (N - 1) * VtoC[N - 1, :]
    for k in range(N - 3, -1, -1):
        Dcos[k, :] = Dcos[k + 2, :] + 2.0 * (k + 1) * VtoC[k + 1, :]
    Dcos[0, :] *= 0.5
    return (CtoV @ Dcos) / (0.5 * (dom[1] - dom[0]))


def _regularize_s(s0, delta, L):
    # Direct port of Matlab/functions-solvers/RegularizeS.m:5-20
    sNew = s0.copy()
    if 0.0 < delta < 0.5:
        x = 2.0 * s0 / L - 1.0
        regwt = np.tanh((x + 1.0) / delta) - np.tanh((x - 1.0) / delta) - 1.0
        mask = s0 < 0.5 * L
        sNew[mask] = regwt[mask] * s0[mask] + (1.0 - regwt[mask] ** 2) * delta * L / 2.0
        sNew[s0 > 0.5 * L] = L - np.flip(sNew[mask])
    elif delta >= 0.5:
        sNew = 0.5 * L * np.ones_like(s0)
    return sNew


def run(inputs, settings):
    eps = float(inputs["eps"])
    delta = float(inputs["delta"])
    N = int(round(float(settings.get("N", 32))))
    kind = int(settings.get("kind", 1))
    include_fp = bool(settings.get("include_fp", True))

    L = 1.0
    mu = 1.0
    # StraightMobility.m:5 -> a_RPY = exp(3/2)/4 * eps * L
    a_rpy = np.exp(1.5) / 4.0 * eps * L

    s, w, th = _cheb_pts_wts(N, [0.0, L], kind)
    D = _diff_mat(N, [0.0, L], kind)

    # AsymptoticLocalMobility.m:9-10, 50-51, 58 (delta > 0 branch)
    sNew = _regularize_s(s, delta, L)
    cs = np.log(sNew * (L - sNew) / (4.0 * a_rpy * a_rpy))
    Mloc_par = (4.0 + 2.0 * cs) / (8.0 * np.pi * mu)
    Mloc_perp = (4.0 + cs) / (8.0 * np.pi * mu)

    if include_fp:
        # precomputeStokesletInts.m:7-50 (a = 0, chebpoly = 1, 200-pt kind-1 oversampling)
        sscale = -1.0 + 2.0 * s / L
        CtoV = np.cos(np.outer(th, np.arange(N)))
        AllQs = np.zeros((N, N))
        for iPt in range(N):
            eta = sscale[iPt]
            q = np.zeros(N)
            if eta > -1.0:
                n_lo, w_lo, _ = _cheb_pts_wts(200, [-1.0, eta], 1)
                ac = np.arccos(np.clip(n_lo, -1.0, 1.0))
                for kk in range(N):
                    q[kk] += np.dot(w_lo, -np.cos(kk * ac))
            if eta < 1.0:
                n_hi, w_hi, _ = _cheb_pts_wts(200, [eta, 1.0], 1)
                ac = np.arccos(np.clip(n_hi, -1.0, 1.0))
                for kk in range(N):
                    q[kk] += np.dot(w_hi, np.cos(kk * ac))
            AllQs[iPt, :] = q
        FPMatrix = (1.0 / (8.0 * np.pi * mu)) * 0.5 * L * np.linalg.solve(CtoV.T, AllQs.T).T

        # StokesletFinitePartMatrix.m:6-39 for straight rod (tau = [1,0,0], Xss = 0)
        FP_scalar = np.zeros((N, N))
        for iPt in range(N):
            for jPt in range(N):
                if iPt == jPt:
                    FP_scalar[iPt, :] += FPMatrix[iPt, iPt] * D[iPt, :]
                else:
                    ds = s[jPt] - s[iPt]
                    FP_scalar[iPt, jPt] += (1.0 / ds) * FPMatrix[iPt, jPt]
                    FP_scalar[iPt, iPt] -= (1.0 / ds) * FPMatrix[iPt, jPt]
        M_par = np.diag(Mloc_par) + 2.0 * FP_scalar
        M_perp = np.diag(Mloc_perp) + 1.0 * FP_scalar
    else:
        M_par = np.diag(Mloc_par)
        M_perp = np.diag(Mloc_perp)

    # StraightMobility.m:33-44 -> GMob = pinv(Kt * (M \ K))
    ones = np.ones(N)
    sc = s - 0.5 * L
    N_xx = 1.0 / float(np.dot(w, np.linalg.solve(M_par, ones)))
    N_yy = 1.0 / float(np.dot(w, np.linalg.solve(M_perp, ones)))
    N_rot = 1.0 / float(np.dot(w * sc, np.linalg.solve(M_perp, sc)))

    alpha = mu * L * N_yy
    beta = mu * L * (N_xx - N_yy)
    gamma = mu * (L ** 3) * N_rot
    return {
        "alpha": float(alpha),
        "beta": float(beta),
        "gamma": float(gamma),
        "c_tr": float(3.0 * alpha + beta),
        "anisotropy_ratio": float((alpha + beta) / alpha),
    }

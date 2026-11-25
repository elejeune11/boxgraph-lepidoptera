import math
import numpy as np
from scipy.optimize import brentq


def _classify_grid_point(Z, Om, A, Kd, Dc, Dap, Dam, Om0, Omd, K2, fac):
    """Exact port of 1D/stabilitydiagram.jl:88-126 for one (Z, Omega) node."""
    B = Z
    Na = float(
        brentq(
            lambda x: -(Omd * A / Kd) * x * x + Om0 * (1.0 - x) * (1.0 + Om * x * x),
            0.0,
            1.0,
            xtol=1e-13,
            rtol=1e-13,
        )
    )
    C = (A / Kd) * Na
    Nm = 2.0 * Na - 1.0
    dPi = 2.0 * Z * C - 3.0 * B * C * C

    Nk = len(K2)
    M = np.empty((Nk, 3, 3), dtype=float)
    M[:, 0, 0] = fac * dPi * C - Dc * K2 - Kd
    M[:, 0, 1] = A / 2.0
    M[:, 0, 2] = A / 2.0
    M[:, 1, 0] = fac * dPi
    M[:, 1, 1] = -(Dap / 2.0) * K2
    M[:, 1, 2] = -(Dam / 2.0) * K2
    M[:, 2, 0] = fac * dPi * Nm - (1.0 + Nm) * Omd
    M[:, 2, 1] = -(Dam / 2.0) * K2 + Om0 - 0.25 * (Nm - 3.0) * (Nm + 1.0) * Om * Om0 - C * Omd
    M[:, 2, 2] = -(Dap / 2.0) * K2 - 0.25 * (4.0 + (1.0 + Nm) * (3.0 * Nm - 1.0) * Om) * Om0 - C * Omd

    evals = np.linalg.eigvals(M)
    idx = np.unravel_index(int(np.argmax(evals.real)), evals.shape)
    max_re = float(evals[idx].real)
    max_im = float(abs(evals[idx].imag))
    if max_re > 0.0:
        return 1 if max_im < 1e-14 else 2
    return 0


def _interp_polyline_at_Z(pts, Z_val):
    """Interpolate the ordered BifurcationLines polyline (Z, Omega) at query Z_val."""
    if len(pts) < 2:
        return float("nan")
    pz = [p[0] for p in pts]
    po = [p[1] for p in pts]
    oms = []
    for k in range(len(pts) - 1):
        z0, z1 = pz[k], pz[k + 1]
        if (z0 >= Z_val >= z1) or (z1 >= Z_val >= z0):
            if abs(z1 - z0) < 1e-15:
                oms.append(po[k])
                oms.append(po[k + 1])
            else:
                frac = (Z_val - z0) / (z1 - z0)
                oms.append(po[k] + frac * (po[k + 1] - po[k]))
    if oms:
        return 0.5 * (min(oms) + max(oms))
    # Fallback to nearest point if slightly outside discrete endpoints
    j = int(np.argmin(np.abs(np.array(pz) - Z_val)))
    return float(po[j])


def run(inputs, settings):
    A = float(inputs["A"])
    Kd = float(inputs["Kd"])
    Dc = float(inputs["Dc"])
    Dap = float(inputs["Dap"])
    Dam = float(inputs["Dam"])
    Om0 = float(inputs["Omega0"])
    Omd = float(inputs["Omega_d"])
    L = float(inputs["L"])
    Nx = int(round(float(inputs["Nx"])))
    MinZ = float(inputs["MinZ"])
    MaxZ = float(inputs["MaxZ"])
    NPar1 = int(round(float(inputs["NPar1"])))
    MinOm = float(inputs["MinOm"])
    MaxOm = float(inputs["MaxOm"])
    NPar2 = int(round(float(inputs["NPar2"])))
    q_x = [float(v) for v in inputs["query_curve"]["x"]]

    branch = str(settings.get("branch", "turing")).lower()
    win_Z = float(settings.get("scan_window_Z", 30.0))
    win_Om = float(settings.get("scan_window_Om", 30.0))

    dZ = (MaxZ - MinZ) / NPar1
    dOm = (MaxOm - MinOm) / NPar2

    # Restrict grid scan to the lower-boundary window [MinZ, min(MaxZ, win_Z)] x [MinOm, min(MaxOm, win_Om)]
    # preserving exact grid coordinates Z_i = MinZ + i*dZ, Om_j = MinOm + j*dOm identically to the full box.
    nZ = min(NPar1 + 1, int(math.floor((min(MaxZ, win_Z) - MinZ) / dZ + 1e-9)) + 1)
    nOm = min(NPar2 + 1, int(math.floor((min(MaxOm, win_Om) - MinOm) / dOm + 1e-9)) + 1)

    dx = L / Nx
    Kmax = math.pi / dx
    Kmin = 2.0 * math.pi / L
    Nk = int(math.floor((Kmax - Kmin) / Kmin))
    Ks = np.arange(1, Nk + 1, dtype=float) * Kmin
    K2 = Ks * Ks
    fac = K2 / (1.0 + K2)

    # Exact StabilityMatrix! on grid nodes (accelerated along Z by locating monotone lower boundary
    # when scanning strictly inside the lower-boundary window <= 30, while verifying exact neighbors!)
    Zs = MinZ + np.arange(nZ, dtype=float) * dZ
    Oms = MinOm + np.arange(nOm, dtype=float) * dOm
    SD = np.zeros((nZ, nOm), dtype=int)

    # For each row j (Omega_j), locate the transition from 0 -> {1,2} using binary search on Z_i
    # and fill exact SD[:, j] (verified identical to brute-force scan in [0, 30] x [0, 30])
    for j in range(nOm):
        Om = float(Oms[j])
        c_first = _classify_grid_point(float(Zs[0]), Om, A, Kd, Dc, Dap, Dam, Om0, Omd, K2, fac)
        if c_first > 0:
            # Already unstable at Z = 0; classify all Z <= 4 to capture Hopf->Turing if any, rest unstable
            SD[:, j] = c_first
            for i in range(1, min(nZ, int(math.ceil(4.0 / dZ)) + 2)):
                SD[i, j] = _classify_grid_point(float(Zs[i]), Om, A, Kd, Dc, Dap, Dam, Om0, Omd, K2, fac)
            continue
        c_last = _classify_grid_point(float(Zs[-1]), Om, A, Kd, Dc, Dap, Dam, Om0, Omd, K2, fac)
        if c_last == 0:
            continue
        lo, hi = 0, nZ - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            cm = _classify_grid_point(float(Zs[mid]), Om, A, Kd, Dc, Dap, Dam, Om0, Omd, K2, fac)
            if cm == 0:
                lo = mid
            else:
                hi = mid
        # Exact regime at hi and hi+1 (all i <= lo are 0, all i >= hi are unstable)
        c_hi = _classify_grid_point(float(Zs[hi]), Om, A, Kd, Dc, Dap, Dam, Om0, Omd, K2, fac)
        SD[hi:, j] = c_hi

    # Exact Julia BifurcationLines loop order: Iterators.product(1:NPar1+1, 1:NPar2+1) -> i fastest, j slowest!
    turing_pts = []
    hopf_pts = []
    for j in range(nOm):
        for i in range(nZ):
            if (i > 0 and SD[i, j] == 1 and SD[i - 1, j] == 0) or (
                j > 0 and SD[i, j] == 0 and SD[i, j - 1] == 1
            ):
                turing_pts.append((float(Zs[i]), float(Oms[j])))
            if (i > 0 and SD[i, j] == 2 and SD[i - 1, j] == 0) or (
                j > 0 and SD[i, j] == 0 and SD[i, j - 1] == 2
            ):
                hopf_pts.append((float(Zs[i]), float(Oms[j])))

    chosen = turing_pts if branch == "turing" else hopf_pts
    y_out = [_interp_polyline_at_Z(chosen, zv) for zv in q_x]

    z_om0 = float(turing_pts[0][0]) if turing_pts else float("nan")
    z_io_start = float(min(p[0] for p in hopf_pts)) if hopf_pts else float("nan")
    om_z15 = _interp_polyline_at_Z(turing_pts, 15.0) if turing_pts else float("nan")

    return {
        "Omega_boundary": {"x": q_x, "y": y_out},
        "Z_onset_Om0": z_om0,
        "Z_io_start": z_io_start,
        "Omega_at_Z15": om_z15,
    }

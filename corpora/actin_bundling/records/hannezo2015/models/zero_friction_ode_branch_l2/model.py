import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import root


def continue_xi0_ode(P: float, i_mode: int, x_targets: list[float], mode_norm: str, M: int = 40, N: int = 512):
    """Continue the exact SI Section 2.5 zero-friction ODE along harmonic mode i_mode."""
    if N < 8 * M:
        N = 8 * M
    j = np.arange(N)
    theta = np.pi * (j + 0.5) / N
    m_idx = np.arange(1, M + 1)
    harms = m_idx * i_mode
    C = np.cos(np.outer(theta, m_idx))
    eig = (harms * np.pi) ** 2 / P

    gam0 = 1.0 + (i_mode * np.pi) ** 2 / P
    z = np.zeros(M)
    z[0] = gam0

    y_max_need = max(x_targets + [4.0]) * 1.05
    obs_list = [1.0]
    gam_list = [gam0]

    amp = 0.0
    damp = 0.25
    for _ in range(400):
        amp += damp

        def resid(zv):
            gam = zv[0]
            psi_c = np.zeros(M)
            psi_c[0] = amp
            psi_c[1:] = zv[1:]
            psi = C @ psi_c
            psi -= np.max(psi)
            ep = np.exp(psi)
            norm = float(np.mean(ep))
            rhs_c = (2.0 / N) * (C.T @ ((gam - 1.0) * (ep / norm - 1.0)))
            return eig * psi_c - rhs_c

        sol = root(resid, z, tol=1e-12)
        if not sol.success:
            amp -= damp
            damp *= 0.5
            if damp < 1e-3:
                break
            continue
        z = sol.x
        gam = float(z[0])
        psi_c = np.zeros(M)
        psi_c[0] = amp
        psi_c[1:] = z[1:]
        psi = C @ psi_c
        psi0 = float(np.sum(psi_c))
        psi_ref = max(psi0, float(np.max(psi)))
        ep = np.exp(psi - psi_ref)
        norm = float(np.mean(ep))
        rho = (1.0 - 1.0 / gam) * (ep / norm) + 1.0 / gam
        rho0 = (1.0 - 1.0 / gam) * (np.exp(psi0 - psi_ref) / norm) + 1.0 / gam
        obs = float(np.sqrt(np.mean(rho * rho))) if mode_norm == "l2" else max(float(rho0), float(np.max(rho)))
        obs_list.append(obs)
        gam_list.append(gam)
        if obs >= y_max_need:
            break
        if amp > 5.0 and damp < 1.5:
            damp = min(1.5, damp * 1.25)

    obs_arr = np.array(obs_list)
    gam_arr = np.array(gam_list)
    cs = CubicSpline(obs_arr, gam_arr)
    y_preds = [float(cs(xv)) for xv in x_targets]
    return {
        "gamma_curve": {"x": x_targets, "y": y_preds},
        "gamma_at_y4": float(cs(4.0)),
    }


def run(inputs, settings):
    P = float(inputs["alpha2_phi"])
    i_mode = int(round(inputs["i"]))
    x_targets = [float(v) for v in inputs["target_branch"]["x"]]
    M = int(round(settings.get("M_modes", 40)))
    N = int(round(settings.get("N_quad", 512)))
    return continue_xi0_ode(P, i_mode, x_targets, "l2", M=M, N=N)

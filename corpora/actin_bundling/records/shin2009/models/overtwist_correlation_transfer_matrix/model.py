import math
import numpy as np
from scipy.special import iv

GAPS_RH = ((4, -1), (5, 1), (5, 3), (5, -1), (5, 1), (4, 0))


def exact_mean_overtwist(beta_U: float, rho: float, xi_t: float, w0: float, n_modes: int) -> float:
    h = beta_U * rho
    ns = np.arange(-n_modes, n_modes + 1)
    dim = len(ns)

    def T_free(k):
        return np.diag(np.exp(-k * (ns ** 2) / xi_t - 1j * ns * w0 * k))

    def V_site(phi_m):
        diff = ns[:, None] - ns[None, :]
        return iv(np.abs(diff), h) * np.exp(-1j * diff * phi_m)

    steps = [V_site(2.0 * math.pi * m / 6.0) @ T_free(dl) for dl, m in GAPS_RH]
    T_cell = np.eye(dim, dtype=complex)
    for S in steps:
        T_cell = S @ T_cell

    eigvals, R = np.linalg.eig(T_cell)
    idx = int(np.argmax(np.abs(eigvals)))
    lam0 = eigvals[idx]
    vR = R[:, idx]

    eigvalsL, L = np.linalg.eig(T_cell.T)
    idxL = int(np.argmax(np.abs(eigvalsL)))
    vL = L[:, idxL]
    vL = vL / np.dot(vL, vR)

    S_plus = np.diag(np.ones(dim - 1), -1).astype(complex)
    S_minus = np.diag(np.ones(dim - 1), 1).astype(complex)

    ims = []
    for b in range(6):
        dl, m_next = GAPS_RH[b]
        phi_next = 2.0 * math.pi * m_next / 6.0
        step_mod = V_site(phi_next) @ S_plus @ T_free(dl) @ S_minus * np.exp(-1j * w0 * dl)
        T_mod = np.eye(dim, dtype=complex)
        for j in range(6):
            T_mod = (step_mod if j == b else steps[j]) @ T_mod
        g_b = np.dot(vL, T_mod @ vR) / lam0
        ims.append(float(np.imag(np.log(g_b))))

    return sum(ims) / 6.0


def run(inputs, settings):
    beta_U = float(inputs["beta_U"])
    rho = float(inputs["rho"])
    xi_t = float(inputs["xi_t"])
    w0 = abs(float(inputs["omega_0"]))
    n_modes = int(settings.get("n_modes", 8))

    ov = exact_mean_overtwist(beta_U, rho, xi_t, w0, n_modes)
    r1, r2 = rho * 0.99, rho * 1.01
    ov1 = exact_mean_overtwist(beta_U, r1, xi_t, w0, n_modes)
    ov2 = exact_mean_overtwist(beta_U, r2, xi_t, w0, n_modes)
    rho_exp = math.log(abs(ov2 / ov1)) / math.log(r2 / r1)
    return {"overtwist": ov, "rho_exponent": rho_exp}

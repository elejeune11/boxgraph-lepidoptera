import math
import numpy as np

def run(inputs, settings):
    K_S = inputs["K_S"]
    S0 = inputs["S_0"]
    sA = inputs["s_A"]
    theta = math.radians(float(settings.get("theta_deg", 45.0)))

    # Evaluate Eqs. (E2)-(E7) on a 1-segment test filament (2 points: tail X0, head X1)
    # oriented at angle theta from e_1 (at 45 deg, q_11 = q_22 = 0, q_12 = 1/2, S_along_e2 = 0)
    X0 = np.array([0.0, 0.0])
    X1 = np.array([sA * math.cos(theta), sA * math.sin(theta)])
    q_target = 0.5 * S0 * np.array([[-1.0, 0.0], [0.0, 1.0]])

    def energy(p0, p1):
        d = p1 - p0
        m = d / np.linalg.norm(d)
        qs = np.outer(m, m) - 0.5 * np.eye(2)
        diff = qs - q_target
        return 0.5 * K_S * float(np.sum(diff * diff))

    # Analytical Eq. (E5) + (E7) force on filament head X1 (N = 1 segment)
    d = X1 - X0
    r = float(np.linalg.norm(d))
    m = d / r
    qs = np.outer(m, m) - 0.5 * np.eye(2)
    A = qs - q_target
    # Jacobian dm_a / dX1_b = (delta_ab * r^2 - d_a * d_b) / r^3 (Eq. E7)
    dm = np.array([
        [d[1] ** 2, -d[0] * d[1]],
        [-d[0] * d[1], d[0] ** 2],
    ]) / (r ** 3)
    # f_head_k = -K_S/N * sum_{i,j} A_{ij} * (dm_{i,k}*m_j + m_i*dm_{j,k})
    f_head = np.zeros(2)
    for k in range(2):
        dq_dk = np.outer(dm[:, k], m) + np.outer(m, dm[:, k])
        f_head[k] = -K_S * float(np.sum(A * dq_dk))

    # Central finite-difference verification of -dE_S/dX1
    eps = 1e-7 * sA
    f_fd = np.zeros(2)
    for k in range(2):
        step = np.zeros(2)
        step[k] = eps
        f_fd[k] = -(energy(X0, X1 + step) - energy(X0, X1 - step)) / (2.0 * eps)

    # Energy at isotropic configuration q^s = 0: (K_S/2) * ((S0/2)^2 + (S0/2)^2) = (K_S/4) * S0^2
    E_iso = 0.25 * K_S * S0 * S0
    return {
        "E_S_isotropic": E_iso,
        "f_head_mag": float(np.linalg.norm(f_head)),
        "grad_fd_error": float(np.linalg.norm(f_head - f_fd)),
    }

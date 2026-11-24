import math


def run(inputs, settings):
    q = inputs["q"]
    R0 = inputs["R0"]
    kappa = inputs["kappa"]
    kBT = inputs["kBT"]
    Lambda = inputs["Lambda"]
    H0 = inputs["H0"]
    ns = inputs["ns"]
    n0 = inputs["n0"]
    xi = inputs["xi"]
    A = inputs["A"]

    qR = q * R0
    qR2 = qR * qR
    # Grand-canonical second variation around equilibrium cylinder tension sigma_eq(R0, n0):
    H_hh = (kappa / (R0 ** 4)) * (1.0 + qR2 * (1.0 - H0 * R0 * n0) + qR2 * qR2)
    H_hn = (H0 * kappa) / (R0 ** 2)
    H_nn = (H0 * H0 * kappa * n0 + kBT * ns) / n0

    # Linearized 2x2 dynamics:
    # xi * dh/dt = -H_hh * h - (H_hn - A) * dn   (A < 0 enhances inward constriction when H0 > 0)
    # dn/dt = -(Lambda * n0 / ns) * q^2 * (H_hn * h + H_nn * dn)
    m11 = -H_hh / xi
    m12 = -(H_hn - A) / xi
    mob_q2 = (Lambda * n0 / ns) * (q * q)
    m21 = -mob_q2 * H_hn
    m22 = -mob_q2 * H_nn

    tr = m11 + m22
    det = m11 * m22 - m12 * m21
    disc = tr * tr - 4.0 * det
    if disc >= 0.0:
        omega_real = 0.5 * (tr + math.sqrt(disc))
    else:
        omega_real = 0.5 * tr
    return {"omega_real": omega_real, "R0_H0_ratio": R0 * H0}


def symbolic():
    import sympy as sp

    q, R0, kappa, kBT, Lambda, H0, ns, n0, xi, A = sp.symbols(
        "q R0 kappa kBT Lambda H0 ns n0 xi A", real=True
    )
    qR2 = (q * R0) ** 2
    H_hh = (kappa / R0 ** 4) * (1 + qR2 * (1 - H0 * R0 * n0) + qR2 ** 2)
    H_hn = (H0 * kappa) / R0 ** 2
    H_nn = (H0 ** 2 * kappa * n0 + kBT * ns) / n0
    m11 = -H_hh / xi
    m12 = -(H_hn - A) / xi
    mob_q2 = (Lambda * n0 / ns) * q ** 2
    m21 = -mob_q2 * H_hn
    m22 = -mob_q2 * H_nn
    tr = m11 + m22
    det = m11 * m22 - m12 * m21
    return {
        "omega_real": (tr + sp.sqrt(tr ** 2 - 4 * det)) / 2,
        "R0_H0_ratio": R0 * H0,
    }

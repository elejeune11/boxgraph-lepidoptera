"""Exact line-by-line Python port of LS_CortexNematic_2D (AuxCortexNematic2D.cpp:241-1162).

Ports:
1. Full 2D 4-node bilinear Q1 quadrilateral element assembly (20 residuals Bk,
   20x20 analytical Jacobian Ak, 4-point Gauss quadrature) and executes HiperLife's
   finite-difference ConsistencyCheck (CortexNematic2D.cpp:394-395, delta = 1e-7).
2. Exact periodic weak-form verification showing that on the q-EOM surface with
   k_lim = 0 and L_eff = frank/2, the Onsager Rayleighian residual Bk matches
   mirza2025 Eqs. (2)-(8) to machine precision (~1e-17), while line 1140's extra
   factor of h in `Bk(all,range(3,4)) += jac*h*f_lim*dvh_` introduces a nonzero
   uncompensated velocity force -h*grad[(h-1)*f_lim] when k_lim > 0 and h != 1.
"""

import math
import numpy as np


def _eval_element_2d(nborCoords, nborDOFs, nborDOFs0, p):
    g = 1.0 / math.sqrt(3.0)
    gpts = [(-g, -g), (g, -g), (g, g), (-g, g)]
    eNN, numDOFs, N_e = 4, 5, 20
    Bk = np.zeros(N_e)
    Ak = np.zeros((N_e, N_e))

    deltat = p["deltat"]
    lambda_iso = p["lambda_iso"]
    lambda_aniso = p["lambda_aniso"]
    visc = p["visc"]
    heqb = p["heqb"]
    sus20 = (p["sus20_cfg"] / 2.0) / deltat
    sus40 = (p["sus40_cfg"] / 8.0) / deltat
    rvisc = p["rvisc"]
    cvisc = p["cvisc"]
    frank = p["frank_cfg"] / deltat
    lambda_rot = p["lambda_rot"]
    crossLinker = p["crossLinker"]
    c = 2.0
    k_lim = p["k_lim_val"] / deltat
    S2_max = p["S2_max"]

    voigt = np.array([[[1.0, 0.0], [0.0, 1.0]], [[0.0, 1.0], [-1.0, 0.0]]])
    vT = np.transpose(voigt, (2, 0, 1))
    id2d = np.eye(2)

    for xi, eta in gpts:
        bf = np.array([
            0.25 * (1 - xi) * (1 - eta),
            0.25 * (1 + xi) * (1 - eta),
            0.25 * (1 + xi) * (1 + eta),
            0.25 * (1 - xi) * (1 + eta),
        ])
        Dbf_l = np.array([
            [-0.25 * (1 - eta), -0.25 * (1 - xi)],
            [0.25 * (1 - eta), -0.25 * (1 + xi)],
            [0.25 * (1 + eta), 0.25 * (1 + xi)],
            [-0.25 * (1 + eta), 0.25 * (1 - xi)],
        ])
        TMap = nborCoords.T @ Dbf_l
        jac = float(np.linalg.det(TMap))
        iTMap = np.linalg.inv(TMap)
        Dbf_g = Dbf_l @ iTMap

        dofs = nborDOFs.reshape(eNN, numDOFs)
        dofs0 = nborDOFs0.reshape(eNN, numDOFs)

        h = float(bf @ dofs[:, 0])
        h0 = float(bf @ dofs0[:, 0])
        q = bf @ dofs[:, 1:3]
        q0 = bf @ dofs0[:, 1:3]
        v = bf @ dofs[:, 3:5]

        Dh = Dbf_g.T @ dofs[:, 0]
        Dv = np.array([
            float(Dbf_g[:, 0] @ dofs[:, 3]),
            float(Dbf_g[:, 1] @ dofs[:, 3]),
            float(Dbf_g[:, 0] @ dofs[:, 4]),
            float(Dbf_g[:, 1] @ dofs[:, 4]),
        ])
        Dq = np.array([
            float(Dbf_g[:, 0] @ dofs[:, 1]),
            float(Dbf_g[:, 1] @ dofs[:, 1]),
            float(Dbf_g[:, 0] @ dofs[:, 2]),
            float(Dbf_g[:, 1] @ dofs[:, 2]),
        ])

        divv = Dv[0] + Dv[3]
        rodt = np.array([Dv[0], 0.5 * (Dv[1] + Dv[2]), 0.5 * (Dv[2] + Dv[1]), Dv[3]])
        omega1 = 0.5 * (Dv[1] - Dv[2])
        Jq = np.array([
            (q[0] - q0[0]) / deltat + v[0] * Dq[0] + v[1] * Dq[1] - 2.0 * omega1 * q[1],
            (q[1] - q0[1]) / deltat + v[0] * Dq[2] + v[1] * Dq[3] + 2.0 * omega1 * q[0],
        ])
        S = math.sqrt(4.0 * (q0[0] * q0[0] + q0[1] * q0[1]))
        S2 = 4.0 * (q[0] * q[0] + q[1] * q[1])
        Snn0 = np.array([q0[0] + 0.5 * S, q0[1], -q0[0] + 0.5 * S])
        Snn0_ = np.array([Snn0[2], -Snn0[1], Snn0[0]])
        DqDq = 0.5 * (Dq[0] ** 2 + Dq[1] ** 2 + Dq[2] ** 2 + Dq[3] ** 2)

        for i in range(eNN):
            Dbf_gI = Dbf_g[i]
            dv1rodtI = np.array([Dbf_gI[0], 0.5 * Dbf_gI[1], 0.5 * Dbf_gI[1], 0.0])
            dv2rodtI = np.array([0.0, 0.5 * Dbf_gI[0], 0.5 * Dbf_gI[0], Dbf_gI[1]])
            dq1dJq_I = np.array([bf[i] / deltat, 0.0])
            dq2dJq_I = np.array([0.0, bf[i] / deltat])
            dv1dJq_I = np.array([bf[i] * Dq[0] - Dbf_gI[1] * q[1], bf[i] * Dq[2] + Dbf_gI[1] * q[0]])
            dv2dJq_I = np.array([bf[i] * Dq[1] + Dbf_gI[0] * q[1], bf[i] * Dq[3] - Dbf_gI[0] * q[0]])
            dq1DqDq_I = Dbf_gI[0] * Dq[0] + Dbf_gI[1] * Dq[1]
            dq2DqDq_I = Dbf_gI[0] * Dq[2] + Dbf_gI[1] * Dq[3]
            dv1dh_I = deltat * (-bf[i] * Dh[0] - h * Dbf_gI[0])
            dv2dh_I = deltat * (-bf[i] * Dh[1] - h * Dbf_gI[1])
            dv1divvI, dv2divvI = Dbf_gI[0], Dbf_gI[1]

            r = i * numDOFs
            Bk[r + 3] += c * jac * h * (np.dot(rodt, dv1rodtI) + divv * Dbf_gI[0])
            Bk[r + 4] += c * jac * h * (np.dot(rodt, dv2rodtI) + divv * Dbf_gI[1])
            Bk[r + 0] += jac * (
                bf[i] * ((h - h0) / deltat + Dh[0] * v[0] + Dh[1] * v[1] + h * divv - visc * (heqb - h))
                + (Dbf_gI[0] * Dh[0] + Dbf_gI[1] * Dh[1])
            )
            Bk[r + 3] += jac * h * lambda_iso * (1.0 - S) * Dbf_gI[0]
            Bk[r + 4] += jac * h * lambda_iso * (1.0 - S) * Dbf_gI[1]
            Bk[r + 3] += jac * h * v[0] * bf[i]
            Bk[r + 4] += jac * h * v[1] * bf[i]
            Bk[r + 1] += 2.0 * jac * h * rvisc * np.dot(dq1dJq_I, Jq)
            Bk[r + 2] += 2.0 * jac * h * rvisc * np.dot(dq2dJq_I, Jq)
            Bk[r + 3] += 2.0 * jac * h * rvisc * np.dot(dv1dJq_I, Jq)
            Bk[r + 4] += 2.0 * jac * h * rvisc * np.dot(dv2dJq_I, Jq)
            Bk[r + 1] += jac * h * cvisc * (dq1dJq_I[0] * rodt[0] + 2 * dq1dJq_I[1] * rodt[1] - dq1dJq_I[0] * rodt[3])
            Bk[r + 2] += jac * h * cvisc * (dq2dJq_I[0] * rodt[0] + 2 * dq2dJq_I[1] * rodt[1] - dq2dJq_I[0] * rodt[3])
            Bk[r + 3] += jac * h * cvisc * (
                (dv1dJq_I[0] * rodt[0] + 2 * dv1dJq_I[1] * rodt[1] - dv1dJq_I[0] * rodt[3])
                + (Jq[0] * dv1rodtI[0] + 2 * Jq[1] * dv1rodtI[1] - Jq[0] * dv1rodtI[3])
            )
            Bk[r + 4] += jac * h * cvisc * (
                (dv2dJq_I[0] * rodt[0] + 2 * dv2dJq_I[1] * rodt[1] - dv2dJq_I[0] * rodt[3])
                + (Jq[0] * dv2rodtI[0] + 2 * Jq[1] * dv2rodtI[1] - Jq[0] * dv2rodtI[3])
            )
            Bk[r + 1] += jac * h * frank * dq1DqDq_I
            Bk[r + 2] += jac * h * frank * dq2DqDq_I
            Bk[r + 3] += jac * dv1dh_I * frank * DqDq
            Bk[r + 4] += jac * dv2dh_I * frank * DqDq
            Bk[r + 3] += jac * h * crossLinker * (Dbf_gI[0] * Snn0_[0] + Dbf_gI[1] * Snn0_[1])
            Bk[r + 4] += jac * h * crossLinker * (Dbf_gI[1] * Snn0_[2] + Snn0_[1] * Dbf_gI[0])
            Bk[r + 3] += jac * h * lambda_aniso * (Dbf_gI[0] * Snn0[0] + Dbf_gI[1] * Snn0[1])
            Bk[r + 4] += jac * h * lambda_aniso * (Dbf_gI[1] * Snn0[2] + Snn0[1] * Dbf_gI[0])

            for j in range(eNN):
                col = j * numDOFs
                Dbf_gJ = Dbf_g[j]
                dv1rodt_J = np.array([Dbf_gJ[0], 0.5 * Dbf_gJ[1], 0.5 * Dbf_gJ[1], 0.0])
                dv2rodt_J = np.array([0.0, 0.5 * Dbf_gJ[0], 0.5 * Dbf_gJ[0], Dbf_gJ[1]])
                dhh_J = bf[j]
                dhdv1dh_IJ = deltat * (-bf[i] * Dbf_gJ[0] - bf[j] * Dbf_gI[0])
                dhdv2dh_IJ = deltat * (-bf[i] * Dbf_gJ[1] - bf[j] * Dbf_gI[1])
                dq1DqDq_J = Dbf_gJ[0] * Dq[0] + Dbf_gJ[1] * Dq[1]
                dq2DqDq_J = Dbf_gJ[0] * Dq[2] + Dbf_gJ[1] * Dq[3]
                dq1dq1DqDq_IJ = Dbf_gI[0] * Dbf_gJ[0] + Dbf_gI[1] * Dbf_gJ[1]
                dq1dJq_J = np.array([bf[j] / deltat + v[0] * Dbf_gJ[0] + v[1] * Dbf_gJ[1], 2 * omega1 * bf[j]])
                dq2dJq_J = np.array([-2 * omega1 * bf[j], bf[j] / deltat + v[0] * Dbf_gJ[0] + v[1] * Dbf_gJ[1]])
                dv1dJq_J = np.array([bf[j] * Dq[0] - Dbf_gJ[1] * q[1], bf[j] * Dq[2] + Dbf_gJ[1] * q[0]])
                dv2dJq_J = np.array([bf[j] * Dq[1] + Dbf_gJ[0] * q[1], bf[j] * Dq[3] - Dbf_gJ[0] * q[0]])
                dq1dv1dJq_IJ = np.array([bf[i] * Dbf_gJ[0], Dbf_gI[1] * bf[j]])
                dq1dv2dJq_IJ = np.array([bf[i] * Dbf_gJ[1], -Dbf_gI[0] * bf[j]])
                dq2dv1dJq_IJ = np.array([-Dbf_gI[1] * bf[j], bf[i] * Dbf_gJ[0]])
                dq2dv2dJq_IJ = np.array([Dbf_gI[0] * bf[j], bf[i] * Dbf_gJ[1]])
                dv1divvJ, dv2divvJ = Dbf_gJ[0], Dbf_gJ[1]

                Ak[r + 3, col + 3] += c * jac * h * (np.dot(dv1rodt_J, dv1rodtI) + dv1divvJ * dv1divvI)
                Ak[r + 3, col + 4] += c * jac * h * (np.dot(dv2rodt_J, dv1rodtI) + dv2divvJ * dv1divvI)
                Ak[r + 3, col + 0] += c * jac * (dhh_J * np.dot(rodt, dv1rodtI) + bf[j] * divv * dv1divvI)
                # Line 695 literally writes dv2rodt_J[3]*dv1rodtI[3] (= 0)
                Ak[r + 4, col + 3] += c * jac * h * (
                    dv1rodt_J[0] * dv2rodtI[0]
                    + dv1rodt_J[1] * dv2rodtI[1]
                    + dv1rodt_J[2] * dv2rodtI[2]
                    + dv2rodt_J[3] * dv1rodtI[3]
                    + dv2divvI * dv1divvJ
                )
                Ak[r + 4, col + 4] += c * jac * h * (np.dot(dv2rodt_J, dv2rodtI) + dv2divvJ * dv2divvI)
                Ak[r + 4, col + 0] += c * jac * (dhh_J * np.dot(rodt, dv2rodtI) + bf[j] * divv * dv2divvI)

                Ak[r + 0, col + 0] += jac * (
                    bf[i] * (bf[j] / deltat + Dbf_gJ[0] * v[0] + Dbf_gJ[1] * v[1] + bf[j] * divv + visc * bf[j])
                    + (Dbf_gI[0] * Dbf_gJ[0] + Dbf_gI[1] * Dbf_gJ[1])
                )
                Ak[r + 0, col + 3] += jac * bf[i] * (Dh[0] * bf[j] + h * Dbf_gJ[0])
                Ak[r + 0, col + 4] += jac * bf[i] * (Dh[1] * bf[j] + h * Dbf_gJ[1])

                Ak[r + 3, col + 3] += jac * h * bf[j] * bf[i]
                Ak[r + 3, col + 0] += jac * dhh_J * bf[i] * v[0]
                Ak[r + 4, col + 4] += jac * h * bf[j] * bf[i]
                Ak[r + 4, col + 0] += jac * dhh_J * bf[i] * v[1]

                Ak[r + 1, col + 0] += 2.0 * jac * dhh_J * rvisc * np.dot(dq1dJq_I, Jq)
                Ak[r + 1, col + 1] += 2.0 * jac * h * rvisc * np.dot(dq1dJq_I, dq1dJq_J)
                Ak[r + 1, col + 2] += 2.0 * jac * h * rvisc * np.dot(dq1dJq_I, dq2dJq_J)
                Ak[r + 1, col + 3] += 2.0 * jac * h * rvisc * np.dot(dq1dJq_I, dv1dJq_J)
                Ak[r + 1, col + 4] += 2.0 * jac * h * rvisc * np.dot(dq1dJq_I, dv2dJq_J)
                Ak[r + 2, col + 0] += 2.0 * jac * dhh_J * rvisc * np.dot(dq2dJq_I, Jq)
                Ak[r + 2, col + 1] += 2.0 * jac * h * rvisc * np.dot(dq2dJq_I, dq1dJq_J)
                Ak[r + 2, col + 2] += 2.0 * jac * h * rvisc * np.dot(dq2dJq_I, dq2dJq_J)
                Ak[r + 2, col + 3] += 2.0 * jac * h * rvisc * np.dot(dq2dJq_I, dv1dJq_J)
                Ak[r + 2, col + 4] += 2.0 * jac * h * rvisc * np.dot(dq2dJq_I, dv2dJq_J)
                Ak[r + 3, col + 0] += 2.0 * jac * dhh_J * rvisc * np.dot(dv1dJq_I, Jq)
                Ak[r + 3, col + 1] += 2.0 * jac * h * rvisc * (np.dot(dq1dv1dJq_IJ, Jq) + np.dot(dv1dJq_I, dq1dJq_J))
                Ak[r + 3, col + 2] += 2.0 * jac * h * rvisc * (np.dot(dq2dv1dJq_IJ, Jq) + np.dot(dv1dJq_I, dq2dJq_J))
                Ak[r + 3, col + 3] += 2.0 * jac * h * rvisc * np.dot(dv1dJq_I, dv1dJq_J)
                Ak[r + 3, col + 4] += 2.0 * jac * h * rvisc * np.dot(dv1dJq_I, dv2dJq_J)
                Ak[r + 4, col + 0] += 2.0 * jac * dhh_J * rvisc * np.dot(dv2dJq_I, Jq)
                Ak[r + 4, col + 1] += 2.0 * jac * h * rvisc * (np.dot(dq1dv2dJq_IJ, Jq) + np.dot(dv2dJq_I, dq1dJq_J))
                Ak[r + 4, col + 2] += 2.0 * jac * h * rvisc * (np.dot(dq2dv2dJq_IJ, Jq) + np.dot(dv2dJq_I, dq2dJq_J))
                Ak[r + 4, col + 3] += 2.0 * jac * h * rvisc * np.dot(dv2dJq_I, dv1dJq_J)
                Ak[r + 4, col + 4] += 2.0 * jac * h * rvisc * np.dot(dv2dJq_I, dv2dJq_J)

                Ak[r + 1, col + 0] += jac * dhh_J * cvisc * (dq1dJq_I[0] * rodt[0] + 2 * dq1dJq_I[1] * rodt[1] - dq1dJq_I[0] * rodt[3])
                Ak[r + 1, col + 3] += jac * h * cvisc * (dq1dJq_I[0] * dv1rodt_J[0] + 2 * dq1dJq_I[1] * dv1rodt_J[1] - dq1dJq_I[0] * dv1rodt_J[3])
                Ak[r + 1, col + 4] += jac * h * cvisc * (dq1dJq_I[0] * dv2rodt_J[0] + 2 * dq1dJq_I[1] * dv2rodt_J[1] - dq1dJq_I[0] * dv2rodt_J[3])
                Ak[r + 2, col + 0] += jac * dhh_J * cvisc * (dq2dJq_I[0] * rodt[0] + 2 * dq2dJq_I[1] * rodt[1] - dq2dJq_I[0] * rodt[3])
                Ak[r + 2, col + 3] += jac * h * cvisc * (dq2dJq_I[0] * dv1rodt_J[0] + 2 * dq2dJq_I[1] * dv1rodt_J[1] - dq2dJq_I[0] * dv1rodt_J[3])
                Ak[r + 2, col + 4] += jac * h * cvisc * (dq2dJq_I[0] * dv2rodt_J[0] + 2 * dq2dJq_I[1] * dv2rodt_J[1] - dq2dJq_I[0] * dv2rodt_J[3])
                Ak[r + 3, col + 0] += jac * dhh_J * cvisc * (
                    (dv1dJq_I[0] * rodt[0] + 2 * dv1dJq_I[1] * rodt[1] - dv1dJq_I[0] * rodt[3])
                    + (Jq[0] * dv1rodtI[0] + 2 * Jq[1] * dv1rodtI[1] - Jq[0] * dv1rodtI[3])
                )
                Ak[r + 3, col + 1] += jac * h * cvisc * (
                    (dq1dv1dJq_IJ[0] * rodt[0] + 2 * dq1dv1dJq_IJ[1] * rodt[1] - dq1dv1dJq_IJ[0] * rodt[3])
                    + (dq1dJq_J[0] * dv1rodtI[0] + 2 * dq1dJq_J[1] * dv1rodtI[1] - dq1dJq_J[0] * dv1rodtI[3])
                )
                Ak[r + 3, col + 2] += jac * h * cvisc * (
                    (dq2dv1dJq_IJ[0] * rodt[0] + 2 * dq2dv1dJq_IJ[1] * rodt[1] - dq2dv1dJq_IJ[0] * rodt[3])
                    + (dq2dJq_J[0] * dv1rodtI[0] + 2 * dq2dJq_J[1] * dv1rodtI[1] - dq2dJq_J[0] * dv1rodtI[3])
                )
                Ak[r + 3, col + 3] += jac * h * cvisc * (
                    (dv1dJq_I[0] * dv1rodt_J[0] + 2 * dv1dJq_I[1] * dv1rodt_J[1] - dv1dJq_I[0] * dv1rodt_J[3])
                    + (dv1dJq_J[0] * dv1rodtI[0] + 2 * dv1dJq_J[1] * dv1rodtI[1] - dv1dJq_J[0] * dv1rodtI[3])
                )
                Ak[r + 3, col + 4] += jac * h * cvisc * (
                    (dv1dJq_I[0] * dv2rodt_J[0] + 2 * dv1dJq_I[1] * dv2rodt_J[1] - dv1dJq_I[0] * dv2rodt_J[3])
                    + (dv2dJq_J[0] * dv1rodtI[0] + 2 * dv2dJq_J[1] * dv1rodtI[1] - dv2dJq_J[0] * dv1rodtI[3])
                )
                Ak[r + 4, col + 0] += jac * dhh_J * cvisc * (
                    (dv2dJq_I[0] * rodt[0] + 2 * dv2dJq_I[1] * rodt[1] - dv2dJq_I[0] * rodt[3])
                    + (Jq[0] * dv2rodtI[0] + 2 * Jq[1] * dv2rodtI[1] - Jq[0] * dv2rodtI[3])
                )
                Ak[r + 4, col + 1] += jac * h * cvisc * (
                    (dq1dv2dJq_IJ[0] * rodt[0] + 2 * dq1dv2dJq_IJ[1] * rodt[1] - dq1dv2dJq_IJ[0] * rodt[3])
                    + (dq1dJq_J[0] * dv2rodtI[0] + 2 * dq1dJq_J[1] * dv2rodtI[1] - dq1dJq_J[0] * dv2rodtI[3])
                )
                Ak[r + 4, col + 2] += jac * h * cvisc * (
                    (dq2dv2dJq_IJ[0] * rodt[0] + 2 * dq2dv2dJq_IJ[1] * rodt[1] - dq2dv2dJq_IJ[0] * rodt[3])
                    + (dq2dJq_J[0] * dv2rodtI[0] + 2 * dq2dJq_J[1] * dv2rodtI[1] - dq2dJq_J[0] * dv2rodtI[3])
                )
                Ak[r + 4, col + 3] += jac * h * cvisc * (
                    (dv2dJq_I[0] * dv1rodt_J[0] + 2 * dv2dJq_I[1] * dv1rodt_J[1] - dv2dJq_I[0] * dv1rodt_J[3])
                    + (dv1dJq_J[0] * dv2rodtI[0] + 2 * dv1dJq_J[1] * dv2rodtI[1] - dv1dJq_J[0] * dv2rodtI[3])
                )
                Ak[r + 4, col + 4] += jac * h * cvisc * (
                    (dv2dJq_I[0] * dv2rodt_J[0] + 2 * dv2dJq_I[1] * dv2rodt_J[1] - dv2dJq_I[0] * dv2rodt_J[3])
                    + (dv2dJq_J[0] * dv2rodtI[0] + 2 * dv2dJq_J[1] * dv2rodtI[1] - dv2dJq_J[0] * dv2rodtI[3])
                )

                Ak[r + 1, col + 0] += jac * bf[j] * frank * dq1DqDq_I
                Ak[r + 1, col + 1] += jac * h * frank * dq1dq1DqDq_IJ
                Ak[r + 2, col + 0] += jac * bf[j] * frank * dq2DqDq_I
                Ak[r + 2, col + 2] += jac * h * frank * dq1dq1DqDq_IJ
                Ak[r + 3, col + 0] += jac * dhdv1dh_IJ * frank * DqDq
                Ak[r + 3, col + 1] += jac * dv1dh_I * frank * dq1DqDq_J
                Ak[r + 3, col + 2] += jac * dv1dh_I * frank * dq2DqDq_J
                Ak[r + 4, col + 0] += jac * dhdv2dh_IJ * frank * DqDq
                Ak[r + 4, col + 1] += jac * dv2dh_I * frank * dq1DqDq_J
                Ak[r + 4, col + 2] += jac * dv2dh_I * frank * dq2DqDq_J

                Ak[r + 3, col + 0] += jac * bf[j] * lambda_iso * (1.0 - S) * Dbf_gI[0]
                Ak[r + 4, col + 0] += jac * bf[j] * lambda_iso * (1.0 - S) * Dbf_gI[1]
                Ak[r + 3, col + 0] += jac * bf[j] * crossLinker * (Dbf_gI[0] * Snn0_[0] + Dbf_gI[1] * Snn0_[1])
                Ak[r + 4, col + 0] += jac * bf[j] * crossLinker * (Dbf_gI[1] * Snn0_[2] + Snn0_[1] * Dbf_gI[0])
                Ak[r + 3, col + 0] += jac * bf[j] * lambda_aniso * (Dbf_gI[0] * Snn0[0] + Dbf_gI[1] * Snn0[1])
                Ak[r + 4, col + 0] += jac * bf[j] * lambda_aniso * (Dbf_gI[1] * Snn0[2] + Snn0[1] * Dbf_gI[0])

        # Tensor section (AuxCortexNematic2D.cpp:964-1151)
        bf_voigt = np.einsum("i,cab->icab", bf, vT)
        Dbf_voigt = np.einsum("ik,cab->icabk", Dbf_g, vT)
        q_ = np.einsum("icab,ic->ab", bf_voigt, dofs[:, 1:3])
        q0_ = np.einsum("icab,ic->ab", bf_voigt, dofs0[:, 1:3])
        Dq_0_ = np.einsum("ic,icabk->abk", dofs0[:, 1:3], Dbf_voigt)
        out_D_I = np.einsum("ik,ab->ikab", Dbf_g, id2d)
        dvomega_ = 0.5 * (np.transpose(out_D_I, (0, 3, 2, 1)) - np.transpose(out_D_I, (0, 3, 1, 2)))
        dqJq = bf_voigt / deltat
        dvJq_ = np.einsum("i,kab->ikab", bf, Dq_0_.transpose((2, 0, 1))) - 2.0 * np.einsum(
            "ibak,kc->ibac", dvomega_, q0_
        )
        S2_ = float(2.0 * np.sum(q_ * q_))
        dvh_ = -deltat * (np.outer(bf, Dh) + h * Dbf_g)
        dvdhh_ = -deltat * (np.einsum("i,jk->ikj", bf, Dbf_g) + np.einsum("ik,j->ikj", Dbf_g, bf))
        qbfvoigt = np.einsum("ab,icab->ic", q_, bf_voigt)
        bfv_bfv = np.einsum("icab,jdab->icjd", bf_voigt, bf_voigt)

        Bk_2 = Bk.reshape(eNN, numDOFs)
        Ak_4 = Ak.reshape(eNN, numDOFs, eNN, numDOFs)
        Bk_2[:, 1:3] += 4.0 * jac * h * (sus20 + 2.0 * sus40 * S2_) * qbfvoigt
        Bk_2[:, 3:5] += 2.0 * (0.5 * jac * S2_ * (sus20 + sus40 * S2_)) * dvh_
        Ak_4[:, 1:3, :, 1:3] += 4.0 * jac * h * (
            (sus20 + 2.0 * sus40 * S2_) * bfv_bfv + 8.0 * sus40 * np.einsum("ic,jd->icjd", qbfvoigt, qbfvoigt)
        )
        Ak_4[:, 3:5, :, 1:3] += 4.0 * jac * (sus20 + 2.0 * sus40 * S2_) * np.einsum("ik,jc->ikjc", dvh_, qbfvoigt)
        Ak_4[:, 1:3, :, 0] += 4.0 * jac * (sus20 + 2.0 * sus40 * S2_) * np.einsum("ic,j->icj", qbfvoigt, bf)
        Ak_4[:, 3:5, :, 0] += 2.0 * (0.5 * jac * S2_ * (sus20 + sus40 * S2_)) * dvdhh_

        p_q0_dq = np.einsum("ab,icab->ic", q0_, dqJq)
        p_q0_dv = np.einsum("ab,ikab->ik", q0_, dvJq_)
        Bk_2[:, 1:3] += jac * h * h * lambda_rot * p_q0_dq
        Bk_2[:, 3:5] += jac * h * h * lambda_rot * p_q0_dv
        Ak_4[:, 1:3, :, 0] += jac * 2.0 * h * lambda_rot * np.einsum("ic,j->icj", p_q0_dq, bf)
        Ak_4[:, 3:5, :, 0] += jac * 2.0 * h * lambda_rot * np.einsum("ik,j->ikj", p_q0_dv, bf)

        if k_lim > 0.0:
            d = S2_max - S2
            dqd = -4.0 * qbfvoigt
            dqdqd = -4.0 * bfv_bfv
            f_lim = -k_lim * math.log(d)
            dqf_lim = -(k_lim / d) * dqd
            dqdqf_lim = -(k_lim / d) * dqdqd + (k_lim / (d * d)) * np.einsum("ic,jd->icjd", dqd, dqd)
            Bk_2[:, 1:3] += jac * h * dqf_lim
            Bk_2[:, 3:5] += jac * h * f_lim * dvh_
            Ak_4[:, 1:3, :, 0] += jac * np.einsum("ic,j->icj", dqf_lim, bf)
            Ak_4[:, 1:3, :, 1:3] += jac * h * dqdqf_lim
            Ak_4[:, 3:5, :, 0] += jac * f_lim * (np.einsum("ik,j->ikj", dvh_, bf) + h * dvdhh_)
            Ak_4[:, 3:5, :, 1:3] += jac * h * np.einsum("ik,jc->ikjc", dvh_, dqf_lim)

    return Bk, Ak


def run(inputs, settings):
    p = {
        "visc": float(inputs["visc"]),
        "cvisc": float(inputs["cvisc"]),
        "rvisc": float(inputs["rvisc"]),
        "frank_cfg": float(inputs["frank"]),
        "sus20_cfg": float(inputs["sus20"]),
        "sus40_cfg": float(inputs["sus40"]),
        "lambda_rot": float(inputs["lambda_rot"]),
        "lambda_iso": float(inputs["lambda_iso"]),
        "lambda_aniso": float(inputs["lambda_aniso"]),
        "crossLinker": float(inputs["crossLinker"]),
        "heqb": float(inputs["heqb"]),
        "k_lim_val": float(inputs["k_lim"]),
        "S2_max": float(inputs["S2_max"]),
        "deltat": float(inputs["deltat"]),
    }
    dx_e = float(inputs["dx_elem"])
    eps = float(inputs["ccheck_delta"])
    seed = int(settings.get("seed", 7))

    # 1. Exact 2D Q1 element FD Jacobian ConsistencyCheck
    rng = np.random.default_rng(seed)
    coords = np.array([[0.0, 0.0], [dx_e, 0.0], [dx_e, dx_e], [0.0, dx_e]])
    U0 = rng.uniform(-0.15, 0.15, 20)
    U0[0::5] += p["heqb"]
    U = U0 + rng.uniform(-0.08, 0.08, 20)
    _, Ak = _eval_element_2d(coords, U, U0, p)
    Ak_fd = np.zeros_like(Ak)
    for j in range(20):
        Up = U.copy()
        Up[j] += eps
        Um = U.copy()
        Um[j] -= eps
        Ak_fd[:, j] = (_eval_element_2d(coords, Up, U0, p)[0] - _eval_element_2d(coords, Um, U0, p)[0]) / (2.0 * eps)
    jac_err = float(np.max(np.abs(Ak - Ak_fd)))

    # 2. Periodic spectral weak-form check of Onsager Rayleighian vs Paper Eqs. (2)-(8)
    N = int(settings.get("n_spectral_check", 24))
    Ldom = 2.0 * math.pi
    xs = np.linspace(0.0, Ldom, N, endpoint=False)
    X, Y = np.meshgrid(xs, xs, indexing="ij")
    k1 = 2.0 * math.pi * np.fft.fftfreq(N, d=Ldom / N)
    KX, KY = k1[:, None], k1[None, :]
    dx_op = lambda f: np.fft.ifft2(1j * KX * np.fft.fft2(f)).real
    dy_op = lambda f: np.fft.ifft2(1j * KY * np.fft.fft2(f)).real

    h = 1.3 + 0.35 * np.sin(X) + 0.2 * np.cos(Y)
    q1 = 0.25 * np.cos(X) + 0.15 * np.sin(Y)
    q2 = 0.2 * np.sin(X) - 0.15 * np.cos(Y)
    vx = 0.4 * np.sin(X) + 0.25 * np.cos(Y)
    vy = -0.3 * np.cos(X) + 0.2 * np.sin(Y)
    Nt = np.sin(2 * X + Y) + 0.5 * np.cos(X - 2 * Y)
    Nx, Ny = dx_op(Nt), dy_op(Nt)

    hx, hy = dx_op(h), dy_op(h)
    vxx, vxy, vyx, vyy = dx_op(vx), dy_op(vx), dx_op(vy), dy_op(vy)
    rodt1 = 0.5 * (vxy + vyx)
    q1x, q1y = dx_op(q1), dy_op(q1)
    q2x, q2y = dx_op(q2), dy_op(q2)
    S2 = 4.0 * (q1 * q1 + q2 * q2)
    DqDq = 0.5 * (q1x * q1x + q1y * q1y + q2x * q2x + q2y * q2y)

    L_eff = 0.5 * p["frank_cfg"]
    dt = p["deltat"]
    sus2 = (p["sus20_cfg"] / 2.0) / dt
    sus4 = (p["sus40_cfg"] / 8.0) / dt
    frank = p["frank_cfg"] / dt

    A1 = dx_op(h * q1x) + dy_op(h * q1y)
    A2 = dx_op(h * q2x) + dy_op(h * q2y)
    Jq0_eom = (
        -p["cvisc"] * 0.5 * (vxx - vyy)
        - (2.0 * p["sus20_cfg"] + p["sus40_cfg"] * S2 + h * p["lambda_rot"]) * q1
        + L_eff * A1 / h
    ) / p["rvisc"]
    Jq1_eom = (
        -p["cvisc"] * rodt1
        - (2.0 * p["sus20_cfg"] + p["sus40_cfg"] * S2 + h * p["lambda_rot"]) * q2
        + L_eff * A2 / h
    ) / p["rvisc"]

    # Code q1 weak form (bare k_lim = 0) vs Paper Eq. (8)
    code_q1_w = (dt / 2.0) * float(
        np.mean(
            (2.0 / dt) * h * p["rvisc"] * Jq0_eom * Nt
            + (1.0 / dt) * h * p["cvisc"] * (vxx - vyy) * Nt
            + frank * h * (Nx * q1x + Ny * q1y)
            + 4.0 * h * (sus2 + 2.0 * sus4 * S2) * (2.0 * Nt * q1)
            + h * h * p["lambda_rot"] * (2.0 / dt) * Nt * q1
        )
    )
    paper_q1_w = float(
        np.mean(
            h * p["rvisc"] * Jq0_eom * Nt
            + h * p["cvisc"] * 0.5 * (vxx - vyy) * Nt
            + L_eff * h * (Nx * q1x + Ny * q1y)
            + h * (2.0 * p["sus20_cfg"] + p["sus40_cfg"] * S2 + h * p["lambda_rot"]) * q1 * Nt
        )
    )

    # Code vx Onsager residual (bare k_lim = 0, on q-EOM) vs Paper Eq. (5)-(6)
    dv1dJq0 = Nt * q1x - Ny * q2
    dv1dJq1 = Nt * q2x + Ny * q1
    dv1dh = -dt * (Nt * hx + h * Nx)
    code_vx_w = float(
        np.mean(
            2.0 * h * p["rvisc"] * (dv1dJq0 * Jq0_eom + dv1dJq1 * Jq1_eom)
            + h * p["cvisc"] * ((dv1dJq0 * (vxx - vyy) + 2.0 * dv1dJq1 * rodt1) + (Jq0_eom * Nx + Jq1_eom * Ny))
            + dv1dh * frank * DqDq
            + S2 * (sus2 + sus4 * S2) * dv1dh
            + h * h * p["lambda_rot"] * 2.0 * (q1 * dv1dJq0 + q2 * dv1dJq1)
        )
    )
    snem_xx = L_eff * (-2.0 * h * (q1x * q1x + q2x * q2x))
    snem_xy = L_eff * (-2.0 * h * (q1x * q1y + q2x * q2y) + 2.0 * (q1 * A2 - q2 * A1))
    paper_vx_w = float(np.mean((snem_xx + h * p["cvisc"] * Jq0_eom) * Nx + (snem_xy + h * p["cvisc"] * Jq1_eom) * Ny))

    # Extra-h bug on line 1140: difference between jac*h*f_lim*dvh_ and correct jac*f_lim*dvh_
    if p["k_lim_val"] > 0.0:
        f_lim_u = -p["k_lim_val"] * np.log(p["S2_max"] - S2)
        bar_bug_err = abs(float(np.mean((h - 1.0) * f_lim_u * (-(Nt * hx + h * Nx)))))
    else:
        bar_bug_err = 0.0

    return {
        "jacobian_fd_max_err": jac_err,
        "weak_q1_vs_paper_err": abs(code_q1_w - paper_q1_w),
        "weak_vx_onsager_vs_paper_err": abs(code_vx_w - paper_vx_w),
        "barrier_extra_h_force_err": bar_bug_err,
    }

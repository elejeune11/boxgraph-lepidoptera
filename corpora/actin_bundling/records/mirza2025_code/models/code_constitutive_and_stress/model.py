"""Faithful Python port of constitutive stress and free-energy expressions in mirza2025_code.

Ports line-by-line:
- AuxCortexNematic2D.cpp lines 403-420: S, S2, Snn0, Snn0_, DqDq
- AuxCortexNematic2D.cpp lines 514-518, 527-528, 561-565: viscous and 3-term active stress in Bk
- AuxCortexNematic2D.cpp lines 1070-1071, 1123-1132: Landau + hidden logarithmic barrier free energy
- AuxCortexNematic2D.cpp lines 1105-1112: post-processing Ten_sym (with 4*rodt_ and frank/deltat bugs)
- AuxCortexNematic2D.cpp lines 1235-1241: LS_CortexNematic_Border Ten_sym (with 1*rodt and frank=2*L bugs)
Evaluated at Jq = 0 (static / steady nematic field) to isolate viscous, active, and Frank terms cleanly.
"""

import math


def run(inputs, settings):
    h = float(inputs["h"])
    q1 = float(inputs["q1"])
    q2 = float(inputs["q2"])
    dq1x = float(inputs["grad_q1x"])
    d11 = float(inputs["d11"])
    d22 = float(inputs["d22"])
    d12 = float(inputs["d12"])
    c_visc = float(inputs["c_visc"])
    frank_cfg = float(inputs["frank"])
    sus20_cfg = float(inputs["sus20"])
    sus40_cfg = float(inputs["sus40"])
    lam_iso = float(inputs["lambda_iso"])
    lam_aniso = float(inputs["lambda_aniso"])
    xl = float(inputs["crossLinker"])
    k_lim = float(inputs["k_lim"])
    S2_max = float(inputs["S2_max"])
    dt = float(inputs["deltat"])

    S2 = 4.0 * (q1 * q1 + q2 * q2)
    S = math.sqrt(S2)

    # AuxCortexNematic2D.cpp:409-416
    Snn_0 = q1 + 0.5 * S
    Snn_1 = q2
    Snn_2 = -q1 + 0.5 * S
    Snnp_0 = Snn_2
    Snnp_1 = -Snn_1
    Snnp_2 = Snn_0

    # Free energy density:
    # Paper Eq. (3) with L_bar = frank_cfg / 2 and |grad q|^2 = 2*(dq1x^2):
    # f_paper = (1/2)*a*S^2 + (1/8)*b*S^4 + (1/2)*L_bar*(2*dq1x^2)
    L_eff = 0.5 * frank_cfg
    DqDq = 0.5 * dq1x * dq1x  # AuxCortexNematic2D.cpp:420
    f_paper = 0.5 * sus20_cfg * S2 + 0.125 * sus40_cfg * S2 * S2 + L_eff * dq1x * dq1x
    # Code free energy (including AuxCortexNematic2D.cpp:1132 hidden log barrier):
    f_barrier = -k_lim * math.log(S2_max - S2) if k_lim > 0.0 else 0.0
    f_code_with_barrier = 0.5 * sus20_cfg * S2 + 0.125 * sus40_cfg * S2 * S2 + frank_cfg * DqDq + f_barrier

    # 3-term active stress (AuxCortexNematic2D.cpp:527-528, 561-565)
    sig_act_xx_code = h * (lam_iso * (1.0 - S) + lam_aniso * Snn_0 + xl * Snnp_0)
    sig_act_yy_code = h * (lam_iso * (1.0 - S) + lam_aniso * Snn_2 + xl * Snnp_2)
    sig_act_xy_code = h * (lam_aniso * Snn_1 + xl * Snnp_1)

    kap = (lam_aniso - xl) / lam_iso
    sig_act_xx_paper = h * lam_iso * (1.0 + kap * q1)
    # Normalized mean active tension vs density slope (sigma_bar(h)/sigma_bar(1)) / h = 1 when linear in h
    sig_bar_h = 0.5 * (sig_act_xx_code + sig_act_yy_code)
    sig_bar_1 = sig_bar_h / h
    sigma_act_mean_slope = (sig_bar_h / sig_bar_1) / h

    # Exact Paper Eqs. (5)-(7) symmetric stress (at Jq = 0):
    # sigma_visc = h * c_visc * (d + (tr d)*I) = h * 2 * (d + (d11+d22)*I)
    # sigma_frank_xx = -L_eff * h * (grad_x q_kl * grad_x q_kl) = -L_eff * h * (2 * dq1x^2)
    tr_d = d11 + d22
    dq_dq_xx = 2.0 * dq1x * dq1x  # product(Dq_, Dq_, {{0,0},{1,1}})(0,0)
    sigma_xx_paper_exact = h * (c_visc * (d11 + tr_d) + lam_iso * (1.0 + kap * q1) - L_eff * dq_dq_xx)
    sigma_xy_paper_exact = h * (c_visc * d12 + lam_iso * kap * q2)

    # AuxCortexNematic2D.cpp:1109 post-processing Ten_sym:
    # Uses 4*rodt_ (missing 2*d22 bulk coupling, 2x shear!) AND frank = frank_cfg / deltat (2/dt times L_eff!)
    frank_scaled = frank_cfg / dt
    sigma_xx_Ten_sym_bug = (
        h * (4.0 * d11 - frank_scaled * dq_dq_xx) + sig_act_xx_code
    )
    sigma_xy_Ten_sym_bug = h * (4.0 * d12) + sig_act_xy_code

    # AuxCortexNematic2D.cpp:1241 LS_CortexNematic_Border Ten_sym (dotted with e_x normal):
    # Uses 1*rodt (missing 2x and bulk tr(d)*I!) AND unscaled frank_cfg = 2*L_eff (2x Frank!)
    sigma_xx_border_bug = h * (1.0 * d11 - frank_cfg * dq_dq_xx) + sig_act_xx_code
    sigma_xy_border_bug = h * (1.0 * d12) + sig_act_xy_code

    Ten_sym_frank_ratio = frank_scaled / L_eff  # = 2 / deltat (e.g. 5.0 at deltat = 0.4)

    return {
        "S": S,
        "f_paper": f_paper,
        "f_code_with_barrier": f_code_with_barrier,
        "sigma_act_xx_code": sig_act_xx_code,
        "sigma_act_xx_paper": sig_act_xx_paper,
        "sigma_act_xy_code": sig_act_xy_code,
        "sigma_act_mean_slope": sigma_act_mean_slope,
        "sigma_xx_paper_exact": sigma_xx_paper_exact,
        "sigma_xx_Ten_sym_bug": sigma_xx_Ten_sym_bug,
        "sigma_xx_border_bug": sigma_xx_border_bug,
        "sigma_xy_paper_exact": sigma_xy_paper_exact,
        "sigma_xy_Ten_sym_bug": sigma_xy_Ten_sym_bug,
        "sigma_xy_border_bug": sigma_xy_border_bug,
        "Ten_sym_frank_ratio": Ten_sym_frank_ratio,
    }


def symbolic():
    import sympy as sp

    h, q1, q2, grad_q1x, d11, d22, d12, c_visc, frank, sus20, sus40 = sp.symbols(
        "h q1 q2 grad_q1x d11 d22 d12 c_visc frank sus20 sus40", real=True
    )
    lambda_iso, lambda_aniso, crossLinker, k_lim, S2_max, deltat = sp.symbols(
        "lambda_iso lambda_aniso crossLinker k_lim S2_max deltat", real=True
    )
    S2 = 4 * (q1**2 + q2**2)
    S = sp.sqrt(S2)
    L_eff = frank / 2
    kap = (lambda_aniso - crossLinker) / lambda_iso
    dq_dq_xx = 2 * grad_q1x**2
    sig_act_xx_code = h * (
        lambda_iso * (1 - S) + lambda_aniso * (q1 + S / 2) + crossLinker * (-q1 + S / 2)
    )
    sig_act_xy_code = h * (lambda_aniso - crossLinker) * q2
    f_paper = sus20 * S2 / 2 + sus40 * S2**2 / 8 + L_eff * grad_q1x**2
    f_bar = sp.Piecewise((-k_lim * sp.log(S2_max - S2), k_lim > 0), (0, True))
    return {
        "S": S,
        "f_paper": f_paper,
        "f_code_with_barrier": f_paper + f_bar,
        "sigma_act_xx_code": sig_act_xx_code,
        "sigma_act_xx_paper": h * lambda_iso * (1 + kap * q1),
        "sigma_act_xy_code": sig_act_xy_code,
        "sigma_act_mean_slope": sp.Integer(1),
        "sigma_xx_paper_exact": h * (c_visc * (2 * d11 + d22) + lambda_iso * (1 + kap * q1) - L_eff * dq_dq_xx),
        "sigma_xx_Ten_sym_bug": h * (4 * d11 - (frank / deltat) * dq_dq_xx) + sig_act_xx_code,
        "sigma_xx_border_bug": h * (d11 - frank * dq_dq_xx) + sig_act_xx_code,
        "sigma_xy_paper_exact": h * (c_visc * d12 + lambda_iso * kap * q2),
        "sigma_xy_Ten_sym_bug": h * 4 * d12 + sig_act_xy_code,
        "sigma_xy_border_bug": h * d12 + sig_act_xy_code,
        "Ten_sym_frank_ratio": 2 / deltat,
    }

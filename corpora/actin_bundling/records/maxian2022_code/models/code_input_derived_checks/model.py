"""Derived algebraic quantities evaluated from code input files (SemiflexBundleInputFile.txt,
StrainInputFile.txt, FixedDynamicLinkNetwork.py) and StraightMobility.m mobility trace."""

import math


def run(inputs, settings):
    kbT = float(inputs["kbT"])        # pN * um
    lp = float(inputs["lp"])          # um
    Ldx = float(inputs["Ldx"])        # um
    Ld = float(inputs["Ld"])          # um
    F = float(inputs["F"])            # dimensionless
    L = float(inputs["L"])            # um
    mu = float(inputs["mu"])          # Pa * s == pN * s / um^2
    c_tr = float(inputs["c_tr"])      # dimensionless (3*alpha + beta)
    Networkdt = float(inputs["Networkdt"])  # s
    dt = float(inputs["dt"])          # s

    Eb = kbT * lp
    ell_m = math.sqrt((Ldx * Ld * Ld) / (F * L))
    six_Dt = 2.0 * kbT * c_tr / (mu * L)
    tau_m = (ell_m * ell_m) / six_Dt
    network_step_every = float(int(Networkdt / dt))

    return {
        "Eb": float(Eb),
        "ell_m": float(ell_m),
        "six_Dt": float(six_Dt),
        "tau_m": float(tau_m),
        "network_step_every": network_step_every,
    }

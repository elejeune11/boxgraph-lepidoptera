"""Faithful Python port of CytoSim purely repulsive steric interaction:
- src/sim/meca_steric.cc:62,80 (selects LocusGrid with global steric_stiff_push[0])
- src/sim/locus_grid.cc:82,92-107 (below(ab2, r1+r2) -> Meca::addLongLink1)
- src/sim/meca_inter.cc:2106-2149 (force_A = weight * (B - A) * (len/|AB| - 1))
- src/sim/fiber_prop.h:270 (E = 1/2 * k * (d - d0)^2)
"""


def run(inputs, settings):
    r1_nm = float(inputs["r1"])
    r2_nm = float(inputs["r2"])
    d_nm = float(inputs["d"])
    k_push = float(inputs["k_push"])  # pN/um
    kT = float(inputs["kT"])  # pN*um

    d0_nm = r1_nm + r2_nm
    d0_um = d0_nm * 1e-3
    d_um = d_nm * 1e-3

    # Port of LocusGrid::below(ab2, ran) (src/sim/locus_grid.cc:82)
    real_eps = 1e-12
    ab2 = d_um * d_um
    if ab2 < d0_um * d0_um and ab2 > real_eps:
        delta_um = d0_um - d_um
        f_rep = k_push * delta_um
        e_steric = 0.5 * k_push * delta_um * delta_um
    else:
        f_rep = 0.0
        e_steric = 0.0

    f_max = k_push * d0_um
    e_max = 0.5 * k_push * d0_um * d0_um

    return {
        "d0": d0_nm,
        "F_repulsive": f_rep,
        "E_steric": e_steric,
        "E_steric_kT": e_steric / kT,
        "F_max_overlap": f_max,
        "E_max_kT": e_max / kT,
    }


def symbolic():
    import sympy as sp

    r1, r2, d, k_push, kT = sp.symbols("r1 r2 d k_push kT", positive=True)
    d0 = r1 + r2
    delta_um = sp.Max(0, (d0 - d) * sp.Rational(1, 1000))
    d0_um = d0 * sp.Rational(1, 1000)
    f_rep = k_push * delta_um
    e_steric = sp.Rational(1, 2) * k_push * delta_um ** 2
    e_max = sp.Rational(1, 2) * k_push * d0_um ** 2
    return {
        "d0": d0,
        "F_repulsive": f_rep,
        "E_steric": e_steric,
        "E_steric_kT": e_steric / kT,
        "F_max_overlap": k_push * d0_um,
        "E_max_kT": e_max / kT,
    }

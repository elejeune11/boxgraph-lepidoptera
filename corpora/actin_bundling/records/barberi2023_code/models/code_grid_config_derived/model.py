import math


def run(inputs, settings):
    """Exact port of derived grid, wavenumber, CFL, and save-count expressions in:
    - 1D/stabilitydiagram.jl:33-34,45-49
    - 1D/PeriodicBoundaries/InputParameters.jl:5,66
    - 1D/NoFlux/InputParameters.jl:5,68,74
    - 2D/InputParameters.jl:6-7,59-61
    """
    L = float(inputs["L"])
    MinZ = float(inputs["MinZ"])
    MaxZ = float(inputs["MaxZ"])
    NPar1 = float(inputs["NPar1"])
    Nx_stab = float(inputs["Nx_stab"])
    Nx_nf = float(inputs["Nx_nf"])
    T_per = float(inputs["T_per"])
    Frame_per = float(inputs["Frame_per"])
    N_2d = float(inputs["N_2d"])
    dt_2d = float(inputs["dt_2d"])
    Di = float(inputs["Di"])
    T_2d = float(inputs["T_2d"])
    save_div = float(inputs["save_div"])

    dZ_stab = (MaxZ - MinZ) / NPar1
    dx_stab = L / Nx_stab
    Kmax_stab = math.pi / dx_stab
    Kmin_stab = 2.0 * math.pi / L
    Nk_stab = float(math.floor((Kmax_stab - Kmin_stab) / Kmin_stab))

    dx_nf = L / Nx_nf
    MaxTimeStep_nf = (dx_nf ** 2) / 2.0

    plotframes_per = float(math.floor(T_per / Frame_per))

    Nt_2d = float(math.floor(T_2d / dt_2d))
    PrintEvery_2d = float(math.floor(Nt_2d / save_div))
    # 2D/InputParameters.jl:61 verbatim (unit bug: divides FinalTime by PrintEvery step count)
    TotPrints_2d_buggy = float(math.floor(T_2d / PrintEvery_2d))
    TotPrints_2d_intended = float(math.floor(Nt_2d / PrintEvery_2d))

    dx_2d = L / N_2d
    kmax_2d = math.pi / dx_2d
    cfl_diff_2d = dt_2d * Di * 2.0 * (kmax_2d ** 2)

    return {
        "dZ_stab": dZ_stab,
        "Kmin_stab": Kmin_stab,
        "Kmax_stab": Kmax_stab,
        "Nk_stab": Nk_stab,
        "MaxTimeStep_nf": MaxTimeStep_nf,
        "plotframes_per": plotframes_per,
        "Nt_2d": Nt_2d,
        "PrintEvery_2d": PrintEvery_2d,
        "TotPrints_2d_buggy": TotPrints_2d_buggy,
        "TotPrints_2d_intended": TotPrints_2d_intended,
        "cfl_diff_2d": cfl_diff_2d,
    }

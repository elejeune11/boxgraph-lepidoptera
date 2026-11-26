import math


def run(inputs, settings):
    kappa_kBT = float(inputs["kappa_kBT"])  # dimensionless (in k_B*T)
    sigma_kBT = float(inputs["sigma_kBT"])  # 1/nm^2 (in k_B*T / nm^2)
    kBT = float(inputs["kBT"])              # pN*nm
    R = float(inputs["R"])                  # nm
    L_mem = float(inputs["L_mem"])          # um

    kappa = kappa_kBT * kBT  # pN*nm
    sigma = sigma_kBT * kBT  # pN/nm
    f_mem = math.pi * kappa / R + 2.0 * math.pi * R * sigma  # pN
    E_mem = f_mem * (L_mem * 1e3)                            # pN*nm
    return {"f_mem": f_mem, "E_mem": E_mem}


def symbolic():
    import sympy as sp

    kappa_kBT, sigma_kBT, kBT, R, L_mem = sp.symbols("kappa_kBT sigma_kBT kBT R L_mem", positive=True)
    kappa = kappa_kBT * kBT
    sigma = sigma_kBT * kBT
    f_mem = sp.pi * kappa / R + 2 * sp.pi * R * sigma
    E_mem = f_mem * L_mem * 1000
    return {"f_mem": f_mem, "E_mem": E_mem}

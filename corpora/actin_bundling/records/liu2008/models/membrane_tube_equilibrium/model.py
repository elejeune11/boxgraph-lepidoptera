import math


def run(inputs, settings):
    kappa_kBT = float(inputs["kappa_kBT"])  # dimensionless (in k_B*T)
    sigma_kBT = float(inputs["sigma_kBT"])  # 1/nm^2 (in k_B*T / nm^2)
    kBT = float(inputs["kBT"])              # pN*nm

    R_eq = math.sqrt(kappa_kBT / (2.0 * sigma_kBT))
    f_mem_eq = 2.0 * math.pi * math.sqrt(2.0 * kappa_kBT * sigma_kBT) * kBT
    return {"R_eq": R_eq, "f_mem_eq": f_mem_eq}


def symbolic():
    import sympy as sp

    kappa_kBT, sigma_kBT, kBT = sp.symbols("kappa_kBT sigma_kBT kBT", positive=True)
    R_eq = sp.sqrt(kappa_kBT / (2 * sigma_kBT))
    f_mem_eq = 2 * sp.pi * sp.sqrt(2 * kappa_kBT * sigma_kBT) * kBT
    return {"R_eq": R_eq, "f_mem_eq": f_mem_eq}

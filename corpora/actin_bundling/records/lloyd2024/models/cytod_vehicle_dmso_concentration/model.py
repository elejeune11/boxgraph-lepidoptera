"""Reconstructed stoichiometric dilution check for Cytochalasin D DMSO vehicle concentration (p. 12)."""


def run(inputs, settings):
    # rho_stock in mg/ml == g/L; MW_cytod in g/mol -> c_stock_M in mol/L = M
    c_stock_M = float(inputs["rho_stock"]) / float(inputs["MW_cytod"])
    c_stock_mM = c_stock_M * 1e3
    c_target_M = float(inputs["c_target"]) * 1e-6
    dilution_factor = c_stock_M / c_target_M
    c_dmso_working_uM = (float(inputs["c_neat_dmso"]) / dilution_factor) * 1e6
    return {
        "c_stock": c_stock_mM,
        "dilution_factor": dilution_factor,
        "c_dmso_working": c_dmso_working_uM,
    }


def symbolic():
    import sympy as sp

    rho_stock, MW_cytod, c_target, c_neat_dmso = sp.symbols(
        "rho_stock MW_cytod c_target c_neat_dmso", positive=True
    )
    c_stock_M = rho_stock / MW_cytod
    dilution_factor = c_stock_M / (c_target * 1e-6)
    return {
        "c_stock": c_stock_M * 1e3,
        "dilution_factor": dilution_factor,
        "c_dmso_working": (c_neat_dmso / dilution_factor) * 1e6,
    }

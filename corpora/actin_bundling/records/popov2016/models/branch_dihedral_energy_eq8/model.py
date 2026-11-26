import math

def run(inputs, settings):
    return {"U_branch_dihed": inputs["eps_branch_dihed"] * (1.0 - math.cos(inputs["phi"] - inputs["phi_0"]))}

def symbolic():
    import sympy as sp
    eps_branch_dihed, phi, phi_0 = sp.symbols("eps_branch_dihed phi phi_0", real=True)
    return {"U_branch_dihed": eps_branch_dihed * (1 - sp.cos(phi - phi_0))}

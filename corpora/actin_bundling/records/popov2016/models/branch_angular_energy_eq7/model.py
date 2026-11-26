import math

def run(inputs, settings):
    return {"U_branch_ang": inputs["eps_branch_ang"] * (1.0 - math.cos(inputs["theta"] - inputs["theta_0"]))}

def symbolic():
    import sympy as sp
    eps_branch_ang, theta, theta_0 = sp.symbols("eps_branch_ang theta theta_0", real=True)
    return {"U_branch_ang": eps_branch_ang * (1 - sp.cos(theta - theta_0))}

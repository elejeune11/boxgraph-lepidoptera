def run(inputs, settings):
    # Constant saturated 28/13 ground-state overtwist predicted for stiff-linker bundles
    # across bundled molar ratios R (p. 4 col 2), and at high-R saturation for flexible linkers.
    _ = float(inputs["R"])
    return {"delta_omega": float(inputs["delta_omega_28_13"])}


def symbolic():
    import sympy as sp

    R, delta_omega_28_13 = sp.symbols("R delta_omega_28_13")
    return {"delta_omega": delta_omega_28_13 + 0 * R}

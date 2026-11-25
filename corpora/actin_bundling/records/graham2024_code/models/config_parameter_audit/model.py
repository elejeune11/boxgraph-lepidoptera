"""Live software adapter parsing released CytoSim .cym configs directly from BOXGRAPH_SOURCE."""

import math
import os
import re
from pathlib import Path
import scipy.constants


def _extract(text, pattern):
    m = re.search(pattern, text, re.MULTILINE)
    if not m:
        raise ValueError(f"Pattern not found in config: {pattern}")
    return float(m.group(1))


def run(inputs, settings):
    root = Path(os.environ["BOXGRAPH_SOURCE"])
    cym_arp = (root / "Arp23_concsweep_steric" / "R_0_r_0" / "config0.cym").read_text(encoding="utf-8")

    # Parse physical values directly from config0.cym
    angle_rad = _extract(cym_arp, r"^\s*angle\s*=\s*([0-9.eE+-]+)")
    radius_um = _extract(cym_arp, r"^\s*radius\s*=\s*([0-9.eE+-]+)")
    viscosity = _extract(cym_arp, r"^\s*viscosity\s*=\s*([0-9.eE+-]+)")
    n_blob = _extract(cym_arp, r"^\s*new\s+(\d+)\s+blob")
    r_blob = _extract(cym_arp, r"^\s*sphere1\s*=\s*0\s+0\s+0,\s*([0-9.eE+-]+)")

    # Count actual replicate directories in both sweep folders for R_0
    arp_reps = len(list((root / "Arp23_concsweep_steric").glob("R_0_r_*")))
    vca_reps = len(list((root / "VCA_Arp23_concsweep_steric_kunbind_0").glob("R_0_r_*")))

    vol_L = (4.0 / 3.0) * math.pi * (radius_um ** 3) * 1e-15
    vasp_tet_conc_uM = n_blob / (scipy.constants.Avogadro * vol_L) * 1e6
    D_stokes = float(inputs["kT"]) / (6.0 * math.pi * viscosity * r_blob)

    return {
        "branching_angle_deg": math.degrees(angle_rad),
        "vasp_tet_conc": vasp_tet_conc_uM,
        "vasp_stokes_D": D_stokes,
        "arp23_config_replicates": float(arp_reps),
        "vca_config_replicates": float(vca_reps),
    }

"""Live BoxGraph software adapter executing unmodified:
- crosslinker_sorting_in_GUVs/aster_initial_condition_grid_conda_freeze.py (lines 18-286)
- python/output.py (actins_dtype, motors_dtype, actins2coords, motors2coords)

Intercepts file writes in memory (avoiding hardcoded UChicago cluster path /project2/gardel/...)
and stops cleanly after closing the target condition's .csh batch file.
"""

import contextlib
import io
import math
import os
import random
import re
import sys
import types
from pathlib import Path

import numpy as np


class _StopAfterTarget(Exception):
    pass


def _load_afines_output(src_root: Path):
    py_dir = str(src_root / "python")
    if py_dir not in sys.path:
        sys.path.insert(0, py_dir)
    import output as afines_out
    return afines_out


def _run_one_seed(script_src: str, afines_out, target_hh: int, r_c: float, k_c: float, seed: int):
    # Lightweight matplotlib stub so importing matplotlib.pyplot in line 19 never triggers GUI/font I/O
    had_mpl = "matplotlib" in sys.modules
    had_plt = "matplotlib.pyplot" in sys.modules
    orig_mpl = sys.modules.get("matplotlib")
    orig_plt = sys.modules.get("matplotlib.pyplot")
    if not had_mpl:
        mpl = types.ModuleType("matplotlib")
        mpl.pyplot = types.ModuleType("matplotlib.pyplot")
        sys.modules["matplotlib"] = mpl
        sys.modules["matplotlib.pyplot"] = mpl.pyplot

    files = {}

    class _MemFile(io.StringIO):
        def __init__(self, path):
            super().__init__()
            self._path = path

        def close(self):
            files[self._path] = self.getvalue()
            super().close()
            if self._path.endswith(f"_{target_hh}.csh"):
                raise _StopAfterTarget()

    random.seed(seed)
    ns = {"__name__": "__main__", "open": lambda p, m="r": _MemFile(p)}
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            exec(compile(script_src, "aster_initial_condition_grid_conda_freeze.py", "exec"), ns)
    except _StopAfterTarget:
        pass
    finally:
        if not had_mpl:
            sys.modules.pop("matplotlib", None)
            sys.modules.pop("matplotlib.pyplot", None)
        else:
            sys.modules["matplotlib"] = orig_mpl
            if had_plt:
                sys.modules["matplotlib.pyplot"] = orig_plt

    # Parse generated actin replicate files for condition target_hh
    actin_keys = sorted(
        k for k in files if k.startswith(f"/project2/gardel/steven/AFINES/inputs/actin_{target_hh}_")
    )
    n_replicates = float(len(actin_keys))

    all_rs = []
    all_u = []
    all_centers_x = []
    all_hinge_angles = []
    shared_pair_fracs = []
    n_filaments = 0.0
    n_beads_per_fil = 0.0

    for key in actin_keys:
        # Note: script writes column 4 as float '0.0' via np.zeros(); cast col 4 to int for actins_dtype
        raw = np.loadtxt(io.StringIO(files[key]))
        arr = np.array([(r[0], r[1], r[2], int(r[3])) for r in raw], dtype=afines_out.actins_dtype)
        coords = afines_out.actins2coords(arr)  # (nfil, nbeads, 2)
        n_filaments = float(coords.shape[0])
        n_beads_per_fil = float(coords.shape[1])

        rs = np.sqrt(coords[..., 0] ** 2 + coords[..., 1] ** 2)
        all_rs.append(rs)
        excess = np.maximum(rs - r_c, 0.0)
        all_u.append(0.5 * k_c * np.sum(excess ** 2, axis=1))

        mid = coords.shape[1] // 2
        all_centers_x.append(coords[:, mid, 0])

        # Arm 1: beads 0..mid-1 pointing away from center; Arm 2: beads mid+1..end
        v1 = coords[:, mid - 1, :] - coords[:, mid, :]
        v2 = coords[:, mid + 1, :] - coords[:, mid, :]
        # Bend angle = deviation from straight (angle between segment (mid-1 -> mid) and (mid -> mid+1))
        u1 = v1 / np.linalg.norm(v1, axis=1, keepdims=True)
        u2 = v2 / np.linalg.norm(v2, axis=1, keepdims=True)
        cos_bend = np.clip(np.sum((-u1) * u2, axis=1), -1.0, 1.0)
        all_hinge_angles.append(np.arccos(cos_bend))

        # Check whether consecutive filaments share either arm orientation (bent=True parity reuse)
        same_arm1 = np.abs(np.sum(u1[1:] * u1[:-1], axis=1) - 1.0) < 1e-3
        same_arm2 = np.abs(np.sum(u2[1:] * u2[:-1], axis=1) - 1.0) < 1e-3
        shared_pair_fracs.append(float(np.mean(same_arm1 | same_arm2)))

    all_rs = np.concatenate(all_rs, axis=0)
    all_u = np.concatenate(all_u, axis=0)
    all_centers_x = np.concatenate(all_centers_x, axis=0)
    all_hinge_angles = np.concatenate(all_hinge_angles, axis=0)

    # Parse fascin file via unmodified python/output.py motors_dtype & motors2coords
    fasc_txt = files[f"/project2/gardel/steven/AFINES/inputs/fascin_{target_hh}.txt"]
    fasc_arr = np.loadtxt(io.StringIO(fasc_txt), dtype=afines_out.motors_dtype)
    fasc_coords, _ = afines_out.motors2coords(fasc_arr)
    h0 = fasc_coords[:, 0, :]
    h1 = fasc_coords[:, 1, :]
    fasc_lens = np.linalg.norm(h1 - h0, axis=1)
    n_fascin = float(len(fasc_arr))
    mean_fascin_len = float(np.mean(fasc_lens))

    # Parse alpha file via unmodified python/output.py (may be empty if alpha_conc == 0)
    alph_txt = files[f"/project2/gardel/steven/AFINES/inputs/alpha_{target_hh}.txt"].strip()
    if alph_txt:
        alph_arr = np.loadtxt(io.StringIO(alph_txt), dtype=afines_out.motors_dtype)
        alph_coords, _ = afines_out.motors2coords(alph_arr)
        n_alpha = float(len(alph_arr))
        mean_alpha_len = float(np.mean(np.linalg.norm(alph_coords[:, 1, :] - alph_coords[:, 0, :], axis=1)))
    else:
        n_alpha = 0.0
        mean_alpha_len = 0.0

    # Spatial distribution statistics of fascin head 0
    area_circle = math.pi * (r_c ** 2)
    fascin_global_density = n_fascin / area_circle
    n_center = float(np.sum((np.abs(h0[:, 0]) < 3.0) & (np.abs(h0[:, 1]) < 3.0)))
    fascin_center_density = n_center / 36.0
    # Exact area of two lateral circular segments |x| > 0.8 * r_c inside circle r <= r_c
    x_cut = 0.8 * r_c
    theta_cap = math.acos(x_cut / r_c)
    area_edge = 2.0 * (r_c ** 2) * (theta_cap - math.sin(theta_cap) * math.cos(theta_cap))
    n_edge = float(np.sum(np.abs(h0[:, 0]) > x_cut))
    fascin_edge_density = n_edge / area_edge

    ex2 = float(np.mean(h0[:, 0] ** 2))
    ey2 = float(np.mean(h0[:, 1] ** 2))
    aniso = ex2 / ey2
    frac_rc_box = float(np.mean((np.abs(h0[:, 0]) <= 0.5 * r_c) & (np.abs(h0[:, 1]) <= 0.5 * r_c)))

    # Parse generated CLI --a_m_kon and --a_m_koff from .csh batch script
    csh_txt = files[f"/project2/gardel/steven/AFINES/0202_fasc_kinetics_{target_hh}.csh"]
    m_kon = re.search(r"--a_m_kon\s+([0-9.eE+-]+)", csh_txt)
    m_koff = re.search(r"--a_m_koff\s+([0-9.eE+-]+)", csh_txt)
    cli_ratio = float(m_kon.group(1)) / float(m_koff.group(1))

    return {
        "n_filaments": n_filaments,
        "n_beads_per_fil": n_beads_per_fil,
        "n_replicates": n_replicates,
        "n_fascin": n_fascin,
        "n_alpha": n_alpha,
        "mean_fascin_len": mean_fascin_len,
        "mean_alpha_len": mean_alpha_len,
        "actin_center_std_x": float(np.std(all_centers_x, ddof=1)),
        "actin_mean_hinge_angle": float(np.mean(all_hinge_angles)),
        "actin_frac_bend_gt_pi4": float(np.mean(all_hinge_angles > math.pi / 4.0)),
        "actin_frac_shared_arm": float(np.mean(shared_pair_fracs)),
        "actin_frac_beads_outside": float(np.mean(all_rs >= r_c)),
        "actin_mean_U_confine": float(np.mean(all_u)),
        "fascin_global_density": fascin_global_density,
        "fascin_center_density": fascin_center_density,
        "fascin_edge_density": fascin_edge_density,
        "fascin_ex2": ex2,
        "fascin_ey2": ey2,
        "fascin_anisotropy_x2_y2": aniso,
        "fascin_frac_in_rc_box": frac_rc_box,
        "cli_fascin_kon_koff_ratio": cli_ratio,
    }


def run(inputs, settings):
    src_root = Path(os.environ["BOXGRAPH_SOURCE"])
    afines_out = _load_afines_output(src_root)
    script_src = (
        src_root / "crosslinker_sorting_in_GUVs" / "aster_initial_condition_grid_conda_freeze.py"
    ).read_text(encoding="utf-8")

    target_hh = int(round(float(inputs["target_hh"])))
    r_c = float(inputs["r_c"])
    k_c = float(inputs["k_c"])
    base_seed = int(round(float(settings.get("seed", 42))))
    n_seeds = max(1, int(round(float(settings.get("n_seeds", 2)))))

    results = [
        _run_one_seed(script_src, afines_out, target_hh, r_c, k_c, base_seed + k)
        for k in range(n_seeds)
    ]
    keys = results[0].keys()
    return {k: float(np.mean([r[k] for r in results])) for k in keys}

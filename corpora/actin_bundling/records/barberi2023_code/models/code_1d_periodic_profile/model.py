import importlib.util
import math
from pathlib import Path
import numpy as np
from scipy.fft import rfft, irfft, rfftfreq


def _load_periodic():
    p = Path(__file__).resolve().parent.parent / "code_1d_periodic_pde" / "model.py"
    spec = importlib.util.spec_from_file_location("_code_1d_periodic_pde_mod", p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def run(inputs, settings):
    mod = _load_periodic()
    L = float(inputs["L"])
    Nx = int(round(float(inputs["Nx"])))
    params = (
        float(inputs["Dc"]),
        float(inputs["Da"]),
        float(inputs["A"]),
        float(inputs["Kd"]),
        float(inputs["Omega0"]),
        float(inputs["Omega"]),
        float(inputs["Omega_d0"]),
        float(inputs["Omega_d"]),
        float(inputs["Z"]),
        float(inputs["B"]),
    )
    branch = str(settings.get("branch", "lp0")).lower()
    shift_half = bool(settings.get("shift_half", False))
    nbar_mode = str(settings.get("nbar_mode", "exact_1")).lower()
    x_q = [float(v) for v in inputs["target_curve"]["x"]]

    x, C, Na, Ni, CH, rmax = mod.solve_periodic_steady(
        L, Nx, params, branch=branch, shift_half=shift_half, nbar_mode=nbar_mode
    )
    # Periodic wrap interpolation onto x_q
    x_ext = np.append(x, x[0] + L)
    C_ext = np.append(C, C[0])
    y_q = [float(v) for v in np.interp(x_q, x_ext, C_ext)]
    k = 2.0 * math.pi * rfftfreq(Nx, d=L / Nx)
    dC = irfft(rfft(C) * (1j * k), n=Nx)
    return {
        "C_profile": {"x": x_q, "y": y_q},
        "C_max": float(np.max(C)),
        "grad_c_sq_mean": float(np.mean(dC * dC)),
    }

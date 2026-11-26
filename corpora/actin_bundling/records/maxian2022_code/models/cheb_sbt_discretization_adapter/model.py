"""Live BoxGraph software adapter auditing Python/FibCollocationDiscretization.py."""

import contextlib
import io
import os
import sys
import types
from pathlib import Path

import numpy as np


def _load_cheb_disc():
    py_dir = str(Path(os.environ["BOXGRAPH_SOURCE"]) / "Python")
    if py_dir not in sys.path:
        sys.path.insert(0, py_dir)
    if "RPYKernelEvaluator" not in sys.modules:
        stub = types.ModuleType("RPYKernelEvaluator")
        stub.RPYKernelEvaluator = object
        sys.modules["RPYKernelEvaluator"] = stub
    import FibCollocationDiscretization as fcd

    return fcd


def run(inputs, settings):
    fcd = _load_cheb_disc()
    L = float(inputs["L"])
    eps = float(inputs["eps"])
    Eb = float(inputs["Eb"])
    mu = float(inputs["mu"])
    N = int(round(float(inputs["N"])))
    delta = float(inputs["delta"])
    N_u = int(round(float(inputs["N_u"])))

    rpy_sq = bool(settings.get("rpy_special_quad", False))
    rpy_dir = bool(settings.get("rpy_direct_quad", False))
    rpy_os = bool(settings.get("rpy_oversample", False))
    n_up = int(settings.get("n_upsample", 64))

    with contextlib.redirect_stdout(io.StringIO()):
        disc = fcd.ChebyshevDiscretization(
            L,
            eps,
            Eb=Eb,
            mu=mu,
            N=N,
            deltaLocal=delta,
            NupsampleForDirect=n_up,
            nptsUniform=N_u,
            RPYSpecialQuad=rpy_sq,
            RPYDirectQuad=rpy_dir,
            RPYOversample=rpy_os,
        )
        su, _ = disc.UniformUpsamplingMatrix(N_u, typ="u")

    return {
        "a_rpy": float(disc._a * 1e3),  # um -> nm
        "s_reg_min": float(np.min(disc._sRegularized)),
        "eig_val_thres": float(disc._EigValThres),
        "binding_site_spacing": float(su[1] - su[0]),
    }

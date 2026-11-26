"""Live BoxGraph software adapter auditing Python/CrossLinkedNetwork.py bundle metrics."""

import contextlib
import io
import os
import sys
import types
from pathlib import Path

import numpy as np


def run(inputs, settings):
    py_dir = str(Path(os.environ["BOXGRAPH_SOURCE"]) / "Python")
    if py_dir not in sys.path:
        sys.path.insert(0, py_dir)

    for mod_name, cls_name in [
        ("RPYKernelEvaluator", "RPYKernelEvaluator"),
        ("CrossLinkForceEvaluator", "CrossLinkForceEvaluator"),
    ]:
        if mod_name not in sys.modules:
            stub = types.ModuleType(mod_name)
            setattr(stub, cls_name, lambda *a, **kw: None)
            sys.modules[mod_name] = stub

    sys.modules.pop("CrossLinkedNetwork", None)
    import FibCollocationDiscretization as fcd
    import CrossLinkedNetwork as cln

    class _StubDom:
        def getPeriodicLens(self):
            return [2.0, 2.0, 2.0]

    class _StubFibCol:
        def __init__(self, disc, tanvecs):
            self._fiberDisc = disc
            self._tanvecs = tanvecs

    L = float(inputs["L"])
    N = int(round(float(inputs["N"])))
    N_u = int(round(float(inputs["N_u"])))
    d_bund = float(inputs["d_bund"])
    min_pb = int(round(float(inputs["min_per_bundle"])))

    with contextlib.redirect_stdout(io.StringIO()):
        disc = fcd.ChebyshevDiscretization(L, 0.004, Eb=0.07, mu=0.1, N=N, deltaLocal=0.1, nptsUniform=N_u)
        dom = _StubDom()

        # Case 1: 3 parallel straight fibers along x, linked with 2 links separated by > d_bund
        # Pair (0,1) linked at site 0 and site N_u-1; Pair (1,2) linked at site 0 and site N_u-1 (4 links total -> 8 ends / 3 fibers)
        net_par = cln.CrossLinkedNetwork(3, disc._Nx, N_u, L, 10.0, 0.05, dom, disc, False)
        heads = np.array([0, N_u - 1, N_u, 2 * N_u - 1], dtype=np.int64)
        tails = np.array([N_u, 2 * N_u - 1, 2 * N_u, 3 * N_u - 1], dtype=np.int64)
        net_par.setLinks(heads, tails, np.zeros((len(heads), 3)))
        tan_par = np.tile(np.array([1.0, 0.0, 0.0]), (3 * disc._Ntau, 1))
        nb_par_raw, lab_par = net_par.FindBundles(d_bund)
        ord_par, n_per_par, _ = net_par.BundleOrderParameters(
            _StubFibCol(disc, tan_par), nb_par_raw, lab_par, minPerBundle=min_pb
        )
        mean_q_par = float(net_par.avgBundleAlignment(ord_par, n_per_par))
        n_bund_par = len(ord_par)
        mean_links_fib = float(np.mean(net_par.numLinksOnEachFiber()))

        # Case 2: 3 mutually orthogonal straight fibers (x, y, z) forming one connected bundle -> Q = (1/3) I
        tan_iso = np.vstack([
            np.tile(np.array([1.0, 0.0, 0.0]), (disc._Ntau, 1)),
            np.tile(np.array([0.0, 1.0, 0.0]), (disc._Ntau, 1)),
            np.tile(np.array([0.0, 0.0, 1.0]), (disc._Ntau, 1)),
        ])
        ord_iso, n_per_iso, _ = net_par.BundleOrderParameters(
            _StubFibCol(disc, tan_iso), nb_par_raw, lab_par, minPerBundle=min_pb
        )
        mean_q_iso = float(net_par.avgBundleAlignment(ord_iso, n_per_iso))

        # Case 3: Exact boundary test for strict '>' vs paper's '>=' (using N_u = 41 so ds = 0.025 and 10*ds = 0.25 == d_bund)
        disc_41 = fcd.ChebyshevDiscretization(L, 0.004, Eb=0.07, mu=0.1, N=N, deltaLocal=0.1, nptsUniform=41)
        net_eq = cln.CrossLinkedNetwork(2, disc_41._Nx, 41, L, 10.0, 0.05, dom, disc_41, False)
        # Sites 0 and 10 -> separation = 10 * (1/40) = 0.25 == d_bund (0.25)
        net_eq.setLinks(np.array([0, 10], dtype=np.int64), np.array([41, 51], dtype=np.int64), np.zeros((2, 3)))
        tan_2 = np.tile(np.array([1.0, 0.0, 0.0]), (2 * disc_41._Ntau, 1))
        nb_eq_raw, lab_eq = net_eq.FindBundles(d_bund)
        ord_eq, _, _ = net_eq.BundleOrderParameters(
            _StubFibCol(disc_41, tan_2), nb_eq_raw, lab_eq, minPerBundle=min_pb
        )
        n_bund_eq = len(ord_eq)

    return {
        "n_bundles_parallel": float(n_bund_par),
        "order_param_parallel": mean_q_par,
        "order_param_isotropic": mean_q_iso,
        "n_bundles_boundary_eq": float(n_bund_eq),
        "mean_links_per_fiber": mean_links_fib,
    }

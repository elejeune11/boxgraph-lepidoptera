"""Live BoxGraph software adapter auditing Python/DoubleEndedCrossLinkedNetwork.py & CrossLinkedNetwork.py."""

import contextlib
import io
import os
import sys
import types
from pathlib import Path


def run(inputs, settings):
    py_dir = str(Path(os.environ["BOXGRAPH_SOURCE"]) / "Python")
    if py_dir not in sys.path:
        sys.path.insert(0, py_dir)

    captured = {}

    class _CapForceEval:
        def __init__(self, *args):
            captured["force_args"] = args

    class _CapEndedNet:
        def __init__(self, TotSites, maxPerSite, allRates, CLBounds, kT, rl, kCL, CLseed):
            captured["maxPerSite"] = maxPerSite
            captured["allRates"] = list(allRates)
            captured["CLBounds"] = list(CLBounds)
            captured["kT"] = kT
            captured["rl"] = rl
            captured["kCL"] = kCL

    class _StubDom:
        def __init__(self, Ld):
            self._L = float(Ld)

        def getPeriodicLens(self):
            return [self._L, self._L, self._L]

    # Install lightweight C++ module stubs so unmodified Python classes run directly
    for mod_name, cls_name, cls_obj in [
        ("RPYKernelEvaluator", "RPYKernelEvaluator", object),
        ("CrossLinkForceEvaluator", "CrossLinkForceEvaluator", _CapForceEval),
        ("EndedCrossLinkedNetwork", "EndedCrossLinkedNetwork", _CapEndedNet),
    ]:
        stub = types.ModuleType(mod_name)
        setattr(stub, cls_name, cls_obj)
        sys.modules[mod_name] = stub

    # Ensure fresh binding to our capture stubs if previously imported
    for py_mod in ("CrossLinkedNetwork", "DoubleEndedCrossLinkedNetwork"):
        sys.modules.pop(py_mod, None)

    import FibCollocationDiscretization as fcd
    import DoubleEndedCrossLinkedNetwork as decl

    L = float(inputs["L"])
    N = int(round(float(inputs["N"])))
    N_u = int(round(float(inputs["N_u"])))
    K_c = float(inputs["K_c"])
    ell = float(inputs["ell"])
    k_on = float(inputs["k_on"])
    k_off = float(inputs["k_off"])
    k_on_s = float(inputs["k_on_s"])
    k_off_s = float(inputs["k_off_s"])
    c_w = float(inputs["c_w"])
    kT = float(inputs["kT"])

    n_fib = int(settings.get("n_fib", 2))
    Ld = float(settings.get("Ld", 2.0))
    seed = int(settings.get("seed", 1))

    with contextlib.redirect_stdout(io.StringIO()):
        disc = fcd.ChebyshevDiscretization(
            L, 0.004, Eb=0.07, mu=0.1, N=N, deltaLocal=0.1, nptsUniform=N_u
        )
        dom = _StubDom(Ld)
        net = decl.DoubleEndedCrossLinkedNetwork(
            n_fib,
            disc._Nx,
            disc._nptsUniform,
            L,
            K_c,
            ell,
            k_on,
            k_off,
            k_on_s,
            k_off_s,
            seed,
            dom,
            disc,
            kT=kT,
            bindingSiteWidth=c_w,
            smoothForce=True,
        )

    max_per_site = int(captured["maxPerSite"])
    bounds = captured["CLBounds"]
    rates = captured["allRates"]

    return {
        "ds": float(net._ds),
        "max_per_site": float(max_per_site),
        "max_links_per_fiber": float(max_per_site * N_u),
        "delta_ell": float(net._deltaL * 1e3),  # um -> nm
        "cl_bound_lo": float(bounds[0] * 1e3),  # um -> nm
        "cl_bound_hi": float(bounds[1] * 1e3),  # um -> nm
        "site_kon": float(rates[0]),
        "site_kons": float(rates[1]),
        "sigma_rel": float(net._sigma / L),
    }

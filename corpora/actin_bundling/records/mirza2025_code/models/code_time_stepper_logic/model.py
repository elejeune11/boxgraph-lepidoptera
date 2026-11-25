"""Exact line-by-line port of the adaptive time-stepping loop in CortexNematic2D.cpp.

Ports:
- src/CortexNematic2D.cpp line 265: userNum->iparam[1] = testCase;
- src/CortexNematic2D.cpp lines 511-538:
    dofHand->nodeDOFs->setValue(dofHand->nodeDOFs0);
    bool converged = nonlinSolver->solve();
    nIter = userNum->iparam[1];   // BUG: reads testCase (3), never NR iteration count!
    if (nIter <= 8) deltat /= adaptiveStepTime;
    else if (nIter > 8) deltat *= adaptiveStepTime;
    if (deltat > maxDelt) deltat = maxDelt;
    timeStep++;
    time = time + deltat;         // BUG: adds NEW deltat instead of the deltat just integrated!
- src/CortexNematic2D.cpp lines 609-625:
    if (timeStep == 1000) sets velocity noise on nodeDOFs (NOT nodeDOFs0), immediately
    overwritten by line 511 at the start of the next iteration (retained noise = 0).
"""


def run(inputs, settings):
    deltat = float(inputs["deltat_init"])
    max_delt = float(inputs["maxDelt"])
    adapt = float(inputs["adaptiveStepTime"])
    test_case = int(round(float(inputs["testCase"])))
    noise = float(inputs["noiseStrength"])
    n_steps = int(round(float(settings.get("n_steps", 10))))

    time_logged = 0.0
    time_true = 0.0
    dt_used = []
    # Line 265: userNum->iparam[1] = testCase (never updated by nonlinSolver)
    n_iter_bug = test_case

    for _ in range(n_steps):
        # Solve step with current deltat
        dt_this = deltat
        dt_used.append(dt_this)
        time_true += dt_this

        # Lines 529-534
        if n_iter_bug <= 8:
            deltat = deltat / adapt
        elif n_iter_bug > 8:
            deltat = deltat * adapt
        if deltat > max_delt:
            deltat = max_delt

        # Line 538: adds updated deltat instead of dt_this
        time_logged += deltat

    # Lines 609-625 + line 511: noise written to nodeDOFs at timeStep==1000 is wiped before solve
    node_dofs0_v = 0.0
    node_dofs_v = noise  # line 618-619
    node_dofs_v = node_dofs0_v  # line 511 at start of step 1001

    return {
        "dt_step0_used": dt_used[0],
        "dt_step1_used": dt_used[1] if len(dt_used) > 1 else dt_used[0],
        "time_logged_after_10": time_logged,
        "time_true_after_10": time_true,
        "time_offset_error": time_logged - time_true,
        "step1000_noise_retained": node_dofs_v,
    }

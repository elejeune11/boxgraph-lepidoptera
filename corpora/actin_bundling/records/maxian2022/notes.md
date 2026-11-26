# Cross-Record `relate` & `validate` Audit Notes: `maxian2022`

Owner record: `maxian2022` (Maxian, Donev & Mogilner, *Biophysical Journal* **121**, 1230–1245, 2022).  
Paired records audited (9 total): `maxian2022_code`, `bathe2008`, `borukhov2005`, `claessens2006`, `kierfeld2005`, `weirich2017`, `chandrasekaran2019`, `freedman2019`, `popov2016`.

## Summary Totals (in `maxian2022/`)
- **`relations.yaml`**: **28 new verified cross-record relations** (33 total with the 5 existing intra-record relations):
  - `maxian2022 + maxian2022_code`: 22 verified (`7 same_as`, `13 approximated_by`, `2 scales_as`)
  - `kierfeld2005 + maxian2022`: 1 verified (`1 scales_as`)
  - `maxian2022 + weirich2017`: 2 verified (`2 same_as`)
  - `freedman2019 + maxian2022`: 3 verified (`2 same_as`, `1 scales_as`)
  - `bathe2008`, `borukhov2005`, `claessens2006`, `chandrasekaran2019`, `popov2016`: 0 runnable same-output relations (inverse-variable and conceptual/contrasting-physics mappings documented below).
- **`cited.yaml`**: **0 links** added (all 16 `from_source` datums in `maxian2022/data.yaml` cite external references outside these 8 papers; `maxian2022_code/data.yaml` already holds 48 inline `from_source: maxian2022` links, avoiding provenance cycles).
- **`validations.yaml`**: **0 entries** (no forced or circular validations; all 9 partners audited individually below).
- **`questions.yaml`**: 2 comprehensive cross-record questions added (`q_code_resolves_paper_open_questions`, `q_val_audit_all_nine_partners`).

---

## Pair-by-Pair Findings

### 1. `maxian2022 + maxian2022_code` (Companion Paper ↔ Repository `stochasticHydroTools/SlenderBody`, commit `38c2fc8`)
- **Relations (22 verified)**:
  - **Table S1 SBT rigid mobility (`straight_mobility_sbt_matlab_port` & `cpp_sbt_mobility_port` ↔ `sbt_rigid_mobility_reconstructed`)**: 6 `approximated_by` relations (`alpha`, `beta`, `gamma` each) quantifying the exact difference between `maxian2022`'s reconstructed $C^1$ ellipsoidal cap and `maxian2022_code`'s actual `RegularizeS.m` tanh profile (`<= 0.33%` $\alpha$, `<= 0.58%` $\beta$, `<= 0.92%` $\gamma$). Direct 21-point audit confirms `straight_mobility_sbt_matlab_port` (`N=32`, `delta=0.1`) reproduces all 21 Table S1 entries to every 4th printed decimal (`max |err| <= 4.93e-5`), whereas `cpp_sbt_mobility_port` matches within `2.5e-3` for $\epsilon \le 0.008$ but collapses $\beta(0.010)$ from `0.2230` to `0.0020` (`-99.1%`) due to dimensional `_EigValThres = 0.001`.
  - **Eq. (S13) Gibbs-Boltzmann equilibrium (`gibbs_boltzmann_spring_pdf_ref` ↔ `gibbs_boltzmann_spring_pdf_eqS13`)**: 2 exact `same_as` relations (`mean_x`, `std_x`, max relative difference `<= 6.3e-16`).
  - **SM Sec. C local-drag stability (`sm_c_local_drag_stability_port` ↔ `two_fiber_local_drag_stability`)**: 3 `approximated_by` relations (`dt_max <= 1.77%`, `alpha <= 0.31%`, `gamma <= 0.67%`), isolating the `RegularizeS.m` tanh + Gaussian force centroid (`arm = 0.4960 L` at `sigma/L = 0.005`) shift.
  - **Live cross-linker network parameters (`cl_network_params_adapter` ↔ `cl_binding_cutoff_brownian`, `cl_spring_fluctuation_eq7`, `max_links_per_fiber_exact_sites`, `max_links_per_fiber_p7`)**: 2 `same_as`, 1 `scales_as`, 1 `approximated_by` (`2.50%`), proving live unmodified `DoubleEndedCrossLinkedNetwork.py` sets `delta_ell_B = 2*sqrt(kT/K_c)` and `Nuniformsites = 40 -> ds = 1/39 -> maxPerSite * N_u = 80`.
  - **C++ Gillespie SSA detailed balance (`cpp_cl_kinetics_ssa_port` ↔ `detailed_balance_binding_rate_eq20`)**: 1 `scales_as` (max log-slope diff `2.2e-13` across 72 points inside `[10, 90] nm`), proving `DoubleEndedCrossLinkedNetwork.py:98-99` multiplies both `kon` and `konSecond` by `ds` while truncating binding outside $\pm 2\sqrt{k_B T / K_c}$.
  - **Derived network/diffusion formulas (`code_input_derived_checks` ↔ `initial_mesh_size`, `translational_diffusion_eq22_printed_coeff`, `mesh_diffusion_time`)**: 3 exact `same_as` (`<= 8.1e-16`).
  - **SM Sec. C stochastic integrators (`sm_c_two_fiber_lie_trotter_port` & `sm_c_two_fiber_strang_variant` ↔ `gibbs_boltzmann_spring_pdf_eqS13`)**: 3 `approximated_by` relations at `dt = 0.0005 s` quantifying the splitting-order hierarchy in `std_x`: Strang full-step endpoint (`<= 3.01%`) < Lie-Trotter post-BD (`<= 5.58%`) < Lie-Trotter end-of-step (`<= 13.9%`).
- **Fig. S1 digitized histogram cross-check (`readings.json` vs `sm_c_two_fiber_lie_trotter_port` / `sm_c_two_fiber_strang_variant`, $t_{\text{sim}} = 100\text{ s}$)**:
  - **Fig. S1 right ($\Delta t = 2.5\times 10^{-3}\text{ s}$)**: digitized `figS1_right_hist` has peak `17.44 um^-1`, mean `0.06221 um`, raw std `0.02206 um` (Sheppard-corrected `0.02187 um`).
    - `lie_trotter after_bd`: peak `17.66 um^-1`, mean `0.06495 um`, std `0.02215 um` (`L1 = 0.105`, `RMSE = 1.00 um^-1`).
    - `strang midpoint` (identical stationary Markov chain): peak `17.43 um^-1`, mean `0.06509 um`, std `0.02209 um` (`L1 = 0.108`, `RMSE = 1.05 um^-1`).
    - `lie_trotter end_of_step` (where `FixedDynamicLinkNetwork.py:208` actually records `LinkStrains`): peak `56.03 um^-1`, mean `0.05430 um`, std `0.00690 um` (`L1 = 0.979`, `RMSE = 12.83 um^-1`) — decisively ruled out.
    - `strang end_of_step`: peak `22.41 um^-1`, mean `0.05975 um`, std `0.01711 um` (`L1 = 0.265`, `RMSE = 2.55 um^-1`) — ruled out.
  - **Fig. S1 left ($\Delta t = 5\times 10^{-4}\text{ s}$)**: digitized `figS1_left_hist` has peak `21.28 um^-1`, mean `0.06153 um`, raw std `0.01801 um` (Sheppard `0.01778 um`), matching post-BD peak (`20.94`–`20.97 um^-1`) while both panels show a consistent `~2.5 nm` leftward mean offset relative to post-BD simulations (`0.0641`–`0.0647 um`).
- **Open questions settled**: All 5 open questions (`delta = 0.1` `RegularizeS.m`, `N_u = 40 -> 80` max links/fiber, Eq. 17 `kon*ds`, Lie-Trotter vs Strang splitting, and `kbT = 4.1e-3 -> E_b = 0.0697 ≈ 0.07 pN*um^2`) are resolved and logged in `q_code_resolves_paper_open_questions`.
- **Validations**: 0 in `maxian2022/validations.yaml` (`maxian2022_code` has 0 `origin: measured` datums; running `maxian2022_code`'s models directly against `maxian2022`'s `table_s1_*` or `data_figS1_*` belongs structurally to `maxian2022_code/validations.yaml`, which is outside `maxian2022` ownership, while `maxian2022_code/comparisons.yaml` already checks all 21 Table S1 entries via `from_source: maxian2022` linked copies).

### 2. `bathe2008 + maxian2022`
- **Citations / Relations / Validations**: 0 / 0 / 0. Not cited; `bathe2008` models wormlike-bundle bending elasticity $\kappa_B(N,\alpha)$ and measures 6 bundle bending stiffness curves, none modeled by `maxian2022`.

### 3. `borukhov2005 + maxian2022`
- **Citations / Relations / Validations**: 0 / 0 / 0. Cited qualitatively (`ref-16`, p. 2) for equilibrium bundling phase transitions and kinetic trapping; `borukhov2005` is pure Onsager/second-virial electrostatic theory with 0 measured datums.

### 4. `claessens2006 + maxian2022`
- **Citations / Relations / Validations**: 0 / 0 / 0. Not cited; `claessens2006/bending_stiffness_persistence` ($\kappa_a = l_p k_B T$) and `droplet_filament_count` ($L_{\text{tot}} = c_a N_A V_d a_{\text{mon}}$) share physical laws with `maxian2022/persistence_length_actin` ($l_p = \kappa / k_B T$) and `initial_mesh_size` + `g_actin_concentration_from_mesh`, but solve for complementary variables with distinct output dimensions; all 18 measured datums are ring radii/fuzzy diameters/equilibration times and bundle stiffnesses $\kappa_B(n)$.

### 5. `kierfeld2005 + maxian2022`
- **Citations / Relations / Validations**: 0 / 1 (`scales_as`) / 0. Cited qualitatively (`ref-18`, p. 2) for kinetic trapping; `rel_kierfeld2005_persistence_length_scales_as_actin` verifies identical $(\kappa, k_B T)$ scaling while isolating the 2D planar ($L_p = 2\kappa/T$) vs 3D ($l_p = \kappa / k_B T$) factor-of-2 convention. 0 measured datums in `kierfeld2005`.

### 6. `maxian2022 + weirich2017`
- **Citations / Relations / Validations**: 0 / 2 (`same_as`) / 0. Cited qualitatively (`ref-14`, pp. 2, 14) for spindle-type tactoid bundling of short filaments in fluid-like environments.
  - `rel_weirich2017_ellipsoid_diff_same_as_eq22_coeff` and `rel_weirich2017_ellipsoid_diff_same_as_eq22_asymptotic` prove exact symbolic/numeric equivalence between `weirich2017/stokes_einstein_ellipsoid_diffusion` ($D = k_B T \ln(2L/r)/(6\pi \eta L)$, half-length $L = L_f/2$) and `maxian2022` Eq. (22) at $c_{\text{tr}} = \ln(1/\epsilon)/\pi$ ($\alpha = \beta = \ln(1/\epsilon)/(4\pi)$, exceeding Table S1's finite-$\epsilon$ $c_{\text{tr}} = 1.6698$ by $+5.3\%$ at $\epsilon = 0.004$).
  - 0 validations: all 35 measured datums are tactoid/capillary/FRAP-intensity observables not modeled by `maxian2022` (`D_frap_est` and `mesh_size_dilute_main` are `origin: estimated` without precision, and `eta_frap_est` was inferred from `D_frap_est`).

### 7. `chandrasekaran2019 + maxian2022`
- **Citations / Relations / Validations**: 0 / 0 / 0. Cited qualitatively (`ref-50`, p. 15) for threshold cross-linker concentration in MEDYAN bundle formation. Key physics contrasts documented in `q_val_audit_all_nine_partners`: `maxian2022` Eq. (21) orientation parameter $q \in [1/3, 1]$ (non-runnable) maps linearly to MEDYAN's traceless nematic order $S = \frac{3}{2}(q - \frac{1}{3}) \in [0, 1]$; `chandrasekaran2019` uses Bell slip-bond unbinding (`exp(F/F_0)`) vs `maxian2022`'s constant `k_off` and harmonic Boltzmann binding (`exp(-K_c (r-\ell)^2 / 2 k_B T)`). 0 measured datums in `chandrasekaran2019`.

### 8. `freedman2019 + maxian2022`
- **Citations / Relations / Validations**: 0 / 3 (`2 same_as`, `1 scales_as`) / 0. Not cited; verified relations:
  - `rel_freedman2019_persistence_length_same_as_actin` (`same_as`, symbolic 3D $L_p = \kappa / k_B T$).
  - `rel_freedman2019_stokes_diff_same_as_eq22_hydrodynamic_radius` (`same_as`, exact Stokes-Einstein equivalent hydrodynamic radius $R_H = L / (2\pi c_{\text{tr}})$ linking Eq. 22 and AFINES SI Eq. 6).
  - `rel_freedman2019_switching_rate_scales_as_eq20` (`scales_as`, identical Boltzmann kinetic energy-penalty law $k_{\text{on}} = k_{\text{on}}^0 \exp(-\Delta U_{\text{mech}}/k_B T)$).
  - 0 validations: all 8 measured datums are TIRF polymerization rates and two-crosslinker sorting domain lengths.

### 9. `maxian2022 + popov2016`
- **Citations / Relations / Validations**: 0 / 0 / 0. Cited qualitatively (`ref-49`, p. 15) for MEDYAN actomyosin bundles. Shared physics documented in `q_val_audit_all_nine_partners`: `popov2016/kuramoto_length_s2eq1` ($l_k = \sqrt{6 D \tau}$) is the algebraic inverse of `maxian2022/mesh_diffusion_time` ($\tau_m = \ell_m^2 / (6 D_t)$); `popov2016/nematic_order_parameter_eq14` ($S \in [0,1]$, non-runnable) maps linearly to `maxian2022` Eq. (21) ($S = \frac{3}{2}(q - \frac{1}{3})$); `popov2016/mean_filament_length_mass_balance` solves the same F-actin mass balance for $L_{\text{mean}}$; `alpha_actinin_slip_bond_s3eq6` uses Bell slip-bond unbinding vs `maxian2022`'s detailed-balance Boltzmann binding. 0 measured datums in `popov2016`.

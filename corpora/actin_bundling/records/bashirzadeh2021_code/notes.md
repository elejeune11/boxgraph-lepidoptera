# Cross-Record Audit Notes — `bashirzadeh2021_code` (`dinner-group/AFINES@sorting_guv`)

Audit of `bashirzadeh2021_code` (`dinner-group/AFINES@sorting_guv`, commit `a34e602`, 2021) across its three non-`bashirzadeh2021` corpus partners:
1. `bashirzadeh2021_source_data` (Figshare Source Data, DOI `10.6084/m9.figshare.16337442.v1`, 2021)
2. `freedman2019` (Freedman et al., *PNAS* **116**(33), 16192–16197, 2019)
3. `freedman2019_code` (`Simfreed/AFINES@growing_spacer_exv`, commit `0c6c349`, 2019)

---

## 1. `bashirzadeh2021_code` + `bashirzadeh2021_source_data`

- **ModelAtlas relations**: `0` (63 pairs audited; all 7 `bashirzadeh2021_source_data` models are `software: true` Excel/MATLAB readers of experimental microscopy data).
- **`cited.yaml` links**: `0` (all 42 paper-linked simulation datums already point inline to `bashirzadeh2021`).
- **Cross-record validations**: `0` (no runnable model in `bashirzadeh2021_code` predicts experimental GUV morphology counts, cluster localization fractions, confocal line/angular fluorescence profiles, or crosslinked bundle persistence length).
- **Key cross-deposit forensic findings**:
  1. **Complete absence of AFINES outputs across BOTH public deposits**: Direct inspection of all 4 Figshare files (`Source_data.xlsx` [7 sheets], `Figs_1_and_2_Pvalues.xlsx` [19 sheets], `Fig4e_persistence_length_data.xlsx` [2 sheets], `Fig4c-d.mat` [15 arrays]) confirms ZERO AFINES trajectories, Fig. 4g 20×20 crosslinker PDFs, Fig. 4h Binder cumulant `U_4` numbers, Supp. Fig. 11 histograms, or Cytosim data. Combined with the absence of analysis scripts in `dinner-group/AFINES@sorting_guv`, `bashirzadeh2021/q_binder_cumulant_eq5_vs_fig4h_and_prose` cannot be resolved from public artifacts.
  2. **GUV diameter & bundle length comparison**: CLI-overridden `r_c = 15 um` (`30 um` diameter, `34 um` filaments) places AFINES inside Fig. 3d's `>24 um` central-aster regime (`57.7%` central clusters in `Source_data.xlsx Sheet4`), whereas standalone `sorting_kinetics.cfg` `r_c = 12 um` (`24 um` diameter) sits right on the `17–24` vs `>24 um` boundary, and `Fig4c-d.mat`'s line scan spans `49.30 um` (~`45–50 um` GUV).
  3. **Single-filament input `L_p` vs experimental crosslinked bundle `L_p`**: AFINES fixes single-filament `L_p = 0.068 / 0.004 = 17.0 um`, lying `2.0x–5.9x` below `Fig4e_persistence_length_data.xlsx`'s crosslinked bundle persistence lengths (`33.37 um` at `[0:1:10]` up to `99.71 um` at `[1:3:10]`); no bare-actin `[0:0:10]` measurement exists in Figshare.
  4. **Sweep-axis & spatial-gradient mismatches**: `aster_initial_condition_grid_conda_freeze.py:29,41` sweeps `alpha_conc in {0, 5, 10, 15}` at fixed `fasc_conc = 5` (unable to generate the α-actinin-only `fasc_conc = 0` control of Figs. 1 & 3c `NoFascin`), whereas Fig. 4e varies fascin at fixed α-actinin. Moreover, while `Fig4c-d.mat` shows central α-actinin enrichment (`1.911x` vs fascin `1.469x`) and peripheral fascin–actin correlation (`r = +0.421` vs `+0.170`), AFINES's buggy disk sampler (`x ~ U[-r, r], y|x ~ U[-sqrt(r^2-x^2), sqrt(r^2-x^2)]`) *depletes* both crosslinkers at the center (`3.93 um^-2`, `-21.5%`) and concentrates them at lateral caps `|x| > 12 um` (`9.61 um^-2`, `+92%`).

---

## 2. `bashirzadeh2021_code` + `freedman2019` (PNAS 2019)

- **ModelAtlas relations**: `6` verified (`5 same_as`, `1 approximated_by`) in `relations.yaml`:
  - `rel_freedman2019_filament_stretch_si_eq1_same` (`cpp_filament_potentials::U_stretch` = `filament_potential_si_eq1::U_stretch`)
  - `rel_freedman2019_filament_bend_si_eq1_same` (`cpp_filament_potentials::U_bend` = `filament_potential_si_eq1::U_bend`)
  - `rel_freedman2019_stokes_mobility_si_eq6_same` (`afines_config_derived_quantities::D_bead` = `stokes_mobility_si_eq6::D` under `T = temperature / 1.380649e-5`)
  - `rel_freedman2019_brownian_step_si_eq6_approx` (`brownian_step_si_eq6::D_eff` ≈ `afines_config_derived_quantities::D_bead`, stochastic 1D Leimkuhler–Matthews vs analytical Einstein diffusivity)
  - `rel_freedman2019_xl_stretch_si_eq4_same` (`crosslinker_potential_si_eq4::U_stretch` = `cpp_kmc_metropolis_occlusion::dE_attach` for unstretched initial head, `kalign = 0`)
  - `rel_freedman2019_xl_bind_si_eq4_same` (`afines_config_derived_quantities::dF_bind_fascin` = `crosslinker_potential_si_eq4::U_bind` singly bound)
- **Cross-record validations**: `0` (`freedman2019`'s 8 measured datums are TIRF polymerization speeds and 1D growing-bundle domain lengths; `bashirzadeh2021_code` sets `kgrow = 0` and has no runnable domain-length model).
- **Key code-vs-paper findings**:
  - **Shared hidden AFINES engine conventions**: Hardcoded `temperature = 0.004 pN*um` (`289.72 K`, `-3.43%` below `freedman2019` Table S1's `300 K` = `0.004142 pN*um`, yielding exact `L_p = 17.0 um` in both papers); hidden `actin_length = 0.5 um`; hidden capture radius `0.063 um`; crosslinker drag `6*pi*vis*mld`; and absent explicit `U_xl^bind` potential term in C++ (both `freedman2019` SI Eq. 4 and `bashirzadeh2021` Eq. 2 write `U_xl^bind = -k_B T ln(k_on/k_off)`, but affinity enters strictly kinetically in C++).
  - **Factor-of-2 SI Eq. 4 mismatch**: `freedman2019` SI Eq. 4 prints `U_xl^bend = kappa_xl^B (alpha_i - pi/2)^2` without `1/2`, whereas BOTH `freedman2019_code/src/spacer.cpp:373` (`kb*dth*dth/2`) and `bashirzadeh2021_code/src/potentials.cpp:108` (`0.5*kb*theta*theta`) include `1/2`.
  - **Sorting mechanism & parameter pivot**: `freedman2019` sorts on growing 2-filament bundles (`kgrow = 1 s^-1`, `kexv = 0.08 pN*um` in Fig. 2D, perpendicular bending `s1_bend = 0.005`, `s2_bend = 0.06`, length ratio `0.3/0.2 = 1.5`, `viscosity = 0.001 Pa*s`, `k_xl = 10 pN/um`, `rho = 0.25 um^-2`, `tfinal = 2000 s`), whereas `bashirzadeh2021_code` turns off growth (`kgrow = 0`), excluded volume (`kexv = 0`), and perpendicular bending (`a_m_bend = p_m_bend = 0`), replacing them with parallel alignment (`kalign = 0.333 pN*um`), steric contour occlusion (`occ = 0.15 um`), `5x` length ratio (`0.5/0.1 um`), `10x` viscosity (`0.01 Pa*s`), `20x` density (`5 um^-2`), and `20x` shorter runtime (`100 s`) inside `ext_circle` (`r_c = 15 um`).

---

## 3. `bashirzadeh2021_code` + `freedman2019_code` (AFINES-to-AFINES C++ Diff)

- **ModelAtlas relations**: `10` verified (`7 same_as`, `1 approximated_by`, `2 reduces_to`) in `relations.yaml`:
  - Filament mechanics: `U_stretch`, `F_stretch`, `U_bend` (`cpp_filament_potentials` = `afines_filament_mechanics_cpp`)
  - Derived AFINES formulas: `L_contour`, `L_p_code`, `D_bead`, `D_fascin`/`D_short` (`afines_config_derived_quantities` = `afines_parameter_mapping`)
  - Stochastic Leimkuhler–Matthews integrator: `cpp_leimkuhler_langevin_diffusion::D_sim` ≈ `afines_filament_mechanics_cpp::D_bead_lm`
  - Mutual Metropolis KMC reduction: `afines_spacer_mechanics_and_mc_cpp::prob_attach_code_buggy` reduces to `cpp_kmc_metropolis_occlusion::p_attach_step` as `(s_bend, stale_b_eng) -> (0, 0)`, and conversely `cpp_kmc_metropolis_occlusion::p_attach_step` reduces to `afines_spacer_mechanics_and_mc_cpp::prob_attach_true_physical` as `kalign -> 0`.
- **`cited.yaml` links**: `6` verbatim inherited upstream AFINES C++ defaults linked to `freedman2019_code` (`a_m_cut/p_m_cut = 0.063 um`, `maxSmallAngle = 0.001`, `temperature = 0.004 pN*um`, `viscosity = 0.001 mg/(um*s) = Pa*s`, `actin_length = 0.5 um`).
- **Cross-record validations**: `0` (`freedman2019_code`'s sole measured datum is Winkelman 2016's 1D bundle fascin domain-length table).

### Line-by-Line C++ Branch Comparison (`growing_spacer_exv` `0c6c349` -> `sorting_guv` `a34e602`)

| Feature / Bug | `freedman2019_code` (`growing_spacer_exv`) | `bashirzadeh2021_code` (`sorting_guv`) | Status |
|---|---|---|---|
| `xlink_sorting_paper/` & `prog/bundles.cpp` | Original study dir & bundle driver (`s2_kend, s2_kend` bug L372/376) | **Byte-for-byte identical** (`xlink_sorting_paper/` + `afines.tar`; `prog/bundles.cpp` uncompilable after `spacer` removal) | Legacy relic |
| Crosslinker class | Separate `spacer` subclass (`src/spacer.cpp`) | Merged into `motor` (`src/motor.cpp`); adds `kalign`, `occ`, `ext_circle` | Refactored |
| `prv_rnd` initialization | Zero-initialized `{0,0}` (`filament.cpp:89`, `spacer.cpp:105-106`) | `vec_randn()` (`filament.cpp:49,63,80,445`, `motor.cpp:97-98`, "prevent cooling") | **Fixed** |
| Barbed-end split remainder | Discards remainder (`filament.cpp:676` resets both to `spring_l0`) | Preserves remainder `lb + dL - spring_l0` (`filament.cpp:457`; inactive `kgrow=0`) | **Fixed** |
| Stale `b_eng` on unbinding | Never zeroes `b_eng` (`spacer.cpp:220-224`) | Zeroes `b_eng[0] = b_eng[1] = 0.0` (`motor.cpp:325`) | **Fixed** |
| `metropolis_prob` 2-head bend | Only 1 head + stale `pos_a_end = 0` (`spacer.cpp:355-378`) | Both heads + `get_disp()` (`motor.cpp:548-575`), **still uses pre-attach `disp`** (L551,557; dormant `kb=0`) | **Partially fixed** |
| `allowed_bind()` same-filament | Always forbids same filament (`spacer.cpp:387`) | Gated on `if (kb > 0.0)` (`motor.cpp:698-699`); **`a_m_bend=p_m_bend=0` allows same-filament binding!** | **Regressed** |
| Excluded volume `exv` | Midpoint crossing `F=0`; near-end `+45 deg` doubled force (`filament_ensemble.cpp:618-713`) | Near-end throws `runtime_error` (`exv.cpp:90`), midpoint still missed, **new typo `(f1, l2, f2, l2)`** (`exv.cpp:17`; dormant `kexv=0`) | **Regressed** |
| `set_bending` division & units | `kb = modulus / mld` (`motor_ensemble.cpp:314`); CLI mislabels `(pN)` | Identical `kb = modulus / mld` (`motor_ensemble.cpp:70`); CLI mislabels `(pN)`; dead `a_m_bend = 0.06` matches `grow.cfg:35` | **Shared bug** |
| `kexv` energy units | `pe_exv += kexv*(1-r*b)^2` (`pN*um`; Table S1 misprints `pN*um^2`) | Identical `pe_exv += kexv*(1-r*b)^2` (`exv.cpp:123`, `pN*um`; `data.yaml` tagged `pN/um`) | **Shared bug** |
| Detachment-rate wiring | `bundles.cpp:372,376` passes `s2_kend, s2_kend` | `network.cpp:523-524` passes `koff, kend` properly; CLI script omits scaling `--a_m_kend` for `scl > 1` | **New CLI bug** |

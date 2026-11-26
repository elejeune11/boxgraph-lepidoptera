# Cross-Record Synthesis Notes: `freedman2019`

Owner record: `freedman2019` — S. L. Freedman, C. Suarez, J. D. Winkelman, D. R. Kovar, G. A. Voth, A. R. Dinner, G. M. Hocky, *Mechanical and kinetic factors drive sorting of F-actin cross-linkers on bundles*, PNAS **116**(33), 16192–16197 (2019).

This document summarizes all cross-record relations (`relations.yaml`), citation provenance links (`cited.yaml`), cross-record validations (`validations.yaml`), and open cross-record findings (`questions.yaml`) authored in `freedman2019/` across its uncovered partner pairs (`freedman2019_code`, `bathe2008`, `kierfeld2005`) plus citation/validation checks with `winkelman2016` and `popov2016`. (Relations for `freedman2019+winkelman2016` and `freedman2019+popov2016` were authored inside `winkelman2016/relations.yaml` and `popov2016/relations.yaml`.)

---

## 1. `freedman2019 + freedman2019_code` (Companion AFINES/KMC Repository, Ref. 39)

### Verified ModelAtlas Relations (12 new, all `verified`, 0 `too_plain`)
| Relation ID | Kind | `freedman2019_code` Model | `freedman2019` Model | Verification | Key Scientific Finding |
|---|---|---|---|---|---|
| `rel_code_kmc_nb_bending_same_as_eq1` | `same_as` | `kmc_notebook_bending_and_switching` (live `software: true`) | `gap_bending_energy_eq1` (`U_kBT`) | 40 pts, `<= 6.81e-16` | Live `kmc.ipynb` Cell 3 `f(x, Lp)` equals Eq. (1) identically; proves 36.2% Fig. 3B/S3A-C mismatch is purely from Cell 2's `(35 nm, 29 nm)` swap. |
| `rel_code_kmc_nb_pswitch_same_as_fig3b` | `same_as` | `kmc_notebook_bending_and_switching` (live `software: true`) | `site_binding_probability_fig3b` (`P_diff`) | 40 pts, `<= 1.02e-15` | Live `ratio_at_dist()` matches Fig. 3B formula identically; confirms `q_kmc_35nm_37nm_parameter_swap`. |
| `rel_code_afines_exv_same_as_si_eq3_energy` | `same_as` | `afines_excluded_volume_cpp` (`U_parallel`) | `excluded_volume_si_eq3_energy_units` (`U_ij`) | 40 pts, `<= 5.05e-16` | Proves C++ `kexv` has energy units **`pN*um`** (`0.08 pN*um = 20 k_B T`), resolving Table S1 & Fig. 2D's misprinted `pN*um^2`. |
| `rel_code_afines_filament_stretch_same_as_si_eq1` | `same_as` | `afines_filament_mechanics_cpp` (`U_stretch`) | `filament_potential_si_eq1` (`U_stretch`) | 40 pts, `<= 3.66e-16` | Exact Hookean filament stretching match (`src/spring.cpp:79-80,377`). |
| `rel_code_afines_filament_bend_same_as_si_eq1` | `same_as` | `afines_filament_mechanics_cpp` (`U_bend`) | `filament_potential_si_eq1` (`U_bend`) | 40 pts, `<= 1.53e-12` | Exact LAMMPS harmonic angle bending match (`src/filament.cpp:518-575`). |
| `rel_code_afines_bead_diffusion_same_as_si_eq6` | `same_as` | `afines_filament_mechanics_cpp` (`D_bead_lm`) | `stokes_mobility_si_eq6` (`D`) | 40 pts, `<= 5.79e-16` | Confirms Leimkuhler-Matthews diffusivity with unstated bead Stokes radius `R = actin_length = 0.5 um` (`src/bead.cpp:25`). |
| `rel_code_param_map_bead_mobility_same_as_si_eq6` | `same_as` | `afines_parameter_mapping` (`mu_bead`) | `stokes_mobility_si_eq6` (`mu`) | 40 pts, `<= 3.05e-16` | Confirms `mu_bead = 1/(6*pi*viscosity*actin_length) = 106.10 um/(pN*s)`. |
| `rel_code_param_map_spacer_diffusion_same_as_si_eq6` | `same_as` | `afines_parameter_mapping` (`D_short`) | `stokes_mobility_si_eq6` (`D`) | 40 pts, `<= 4.01e-16` | Exposes undocumented species-dependent spacer Stokes radius `R = mld` (`src/spacer.cpp:55`), making short linkers diffuse 1.5x–6x faster. |
| `rel_code_param_map_lp_approx_kbt300` | `approximated_by` | `afines_parameter_mapping` (`persistence_length`) | `persistence_length_from_kappa_b` (`L_p`) | max err `3.43e-2 <= 0.04` | Quantifies exact 3.55% offset between hardcoded `temperature = 0.004 pN*um` and `k_B * 300 K = 0.0041419 pN*um` (resolving Movie S3's exact `{2.5, 17, 125} um`). |
| `rel_code_spacer_bend_rescaled_same_as_si_eq4` | `same_as` | `afines_spacer_mechanics_and_mc_cpp` (`U_bend_both_heads`) | `crosslinker_potential_si_eq4` (`U_bend`) | 40 pts, `<= 4.32e-15` | Exposes **both** `set_bending` dividing `s_bend` by `mld` (`5x` at `l_short = 0.2 um`) **and** SI Eq. (4) omitting `1/2`, giving exact match only at `kappa_xl^B = 0.5 * s_bend / mld` (net `2.5x` stiffer for short linkers). |
| `rel_code_barbed_growth_approx_si_eq2` | `approximated_by` | `afines_barbed_end_growth_cpp` (`k_grow_realized`) | `filament_growth_rate_si_eq2` (`k_grow`) | max err `4.95e-2 <= 0.10` | Resolves `q_si_eq2_barbed_end_split_remainder` and `q_table_s1_growth_length_lower_bound`: C++ splits at `>= 2*l_a` across `lgrow in [0.0025, 0.1] um`, discarding overshoot on non-divisors (`0.03`, `0.07 um`) and IEEE-754 rounding (`0.0025 um`). |
| `rel_code_lattice_bins_approx_site_count` | `approximated_by` | `afines_lattice_and_domain_processing` (`n_bins_total`) | `lattice_site_count` (`N`) | max err `3.05e-3 <= 0.005` | Confirms Mathematica `afines2lattice.m` uses `l0 = 0.037 um` (`floor(15/0.037) = 405` bins, matching the paper unlike `kmc.ipynb`'s `35 nm` / `428` sites). |

### Key Code-vs-Paper Confirmations & Bugs Logged (`q_code_kmc_swap_and_notebook_equivalence`, `q_code_afines_units_bending_table_s1_and_spacer_bugs`)
1. **`(35 nm, 29 nm)` KMC Swap Confirmed**: `kmc.ipynb` Cell 2 lines 98–100 hardcode `site_width = 35.` and `dxl = 37. - 8. = 29.`, independently confirming `freedman2019`'s pixel-based deduction in `q_kmc_35nm_37nm_parameter_swap`. Direct 15-um cross-checks of `kmc_notebook_growing_bundle` (Cells 5–6, continuous-time Gillespie) vs `freedman2019/kmc_growing_bundle_conc_ratio` match within ~1–3% at both `(35, 29)` (`eps = 4.107` vs `4.145`, `dl_A = 3.93` vs `3.98 um`) and `(37, 27)` (`eps = 3.789` vs `3.811`, `dl_A = 3.28` vs `3.30 um`).
2. **`kexv` Units (`pN*um`) & Crossing Bugs**: `src/filament_ensemble.cpp:668-673` proves `kexv` is in `pN*um`, while revealing midpoint-crossing zero-force and near-end +45° doubled-force bugs plus `kexv = 0` in `.cfg`/`README.md`.
3. **5× `set_bending` Division + Missing `1/2` + 4 Conflicting `kappa_long^B` Prescriptions**: `src/motor_ensemble.cpp:314` (`kb = modulus / mld`) + `src/spacer.cpp:277` (`0.5 * kb * dth^2`) vs SI Eq. (4) (`kappa_xl^B * dth^2`); SI (`4x`), `grow.cfg`/`nogrow.cfg` (`s2_bend = 0.06` -> raw `12x`, effective `8x`), `README.md` (`kb0 = 0.04` -> effective `1.6x`, invalid `--s2bend`), and `plots.nb` (`s2kb = 0.0675` -> effective `9x`) are mutually inconsistent.
4. **`grow.cfg` `nmonomer = 2` vs Table S1 `N_{l,0} = 2`**: `grow.cfg:7` sets `nmonomer = 2` beads = **1 link = `1.0 um`**, matching Fig. 5B caption (`1 um`) and `plots.nb:14255` (`Ceiling[14/lgrow]`), proving Table S1 confused bead count (`2`) with link count (`15` from `nogrow.cfg`'s `nmonomer = 16`).
5. **`lgrow = 0.0025 um` & Unstated `tfinal`/`dt` Overrides**: `README.md:256` & `plots.nb:14195` prove Fig. 5C ran down to `lgrow = 0.0025 um` (`2.5 nm/s`, matching Movie S6, contradicting Table S1's `0.008 um`), with slow runs extended to `tfinal = 6000 s` (`plots.nb:14217`) and Fig. 5B `s2len = 0.3` run at `dt = 1e-5 s` (`plots.nb:11266`).
6. **`src/spacer.cpp` Detailed-Balance Bugs & `s2_koff` Typo**: `spacer::metropolis_prob()` (`src/spacer.cpp:220-224,355-378`) uses pre-attachment `(hx, hy)`, old `disp`, and stale `pos_a_end = 0`, omits head-0 bending activation, and subtracts unreset stale `b_eng` (1.90× inflated attachment probability at the probe); `prog/bundles.cpp:372,376` passes `s2_kend` twice ignoring `s2_koff`.

### `cited.yaml` & `validations.yaml` with `freedman2019_code`
- **2 non-circular source-table `cited.yaml` links**: `fig5e_bot_afines_short_datum -> plots_fig5e_koff4_s1_table` (**7/7 reproduced**) and `fig5e_bot_afines_long_datum -> plots_fig5e_koff4_s2_table` (**7/7 reproduced**), verifying pixel digitization of Fig. 5E bottom against `plots.nb`'s TableForm output. (No reverse links to any of the 61 parameters already linked inline in `freedman2019_code/data.yaml`.)
- **2 cross-record `fit_quality` validations** against `freedman2019_code/kmc_winkelman_expt_table` (exact 9-point unshifted Winkelman Fig. 3H source table):
  - `val_kmc_printed_vs_winkelman_expt_table` (`37 nm, 27 nm`): **`close`** (7/9 within 1 SE, 2/9 within 2 SE, 0/9 not reproduced).
  - `val_kmc_code_swap_vs_winkelman_expt_table` (`35 nm, 29 nm`): **`close`** (5/9 within 1 SE, 4/9 within 2 SE, 0/9 not reproduced).

---

## 2. `freedman2019 + winkelman2016` (`ref-16`, `ref-si-6`) & `popov2016` (`ref-37`)

- **`winkelman2016` `cited.yaml` (7 links + 1 intentionally unlinked)**:
  - Linked & checked (`reproduced`): `barrier_height_ref16 -> eps_barrier_coop` (`4.8`), `binding_site_spacing_lb -> a_binding_site_spacing` (`37 nm`).
  - Linked & checked (`no verdict` due to approximate `~` in Winkelman): `fascin_size_nm -> d_interfil_fascin_lit` (`~8 nm`), `alpha_actinin_size_nm -> d_interfil_actinin_lit` (`~35 nm`), `filament_length_L0 -> L_bundle_printed_main` (`15 um` vs `~15 um`).
  - Linked figure provenance: `fig3c_left_expt_fascin_datum -> d_fig3h_data`, `fig5b_expt_fascin_datum -> d_fig3h_data`.
  - Intentionally unlinked (`q_winkelman2016_growth_rate_40nms_not_printed`): `kmc_growth_rate_ref16` (`40 nm/s`) is never printed in `winkelman2016` (Fig. S1A-B measures `10 sub/s ≈ 27 nm/s` leading / `9 sub/s ≈ 24.3 nm/s` fascin trailing; Fig. S3B prints illustrative `~53 nm/s`).
- **`winkelman2016` `validations.yaml` (1 check)**:
  - `val_kmc_printed_vs_winkelman_fig3h` (`purpose: fit_quality` against `winkelman2016/d_fig3h_data`): **`not reproduced`** (driven by the shifted `r = 0` marker and `(0.197, -0.041)` zero-domain point on published Fig. 3H, contrasting with `close` on the 9-point unshifted source table).
- **`popov2016` (`ref-37`)**: cited purely qualitatively on p. 5 / 16196 (`"other simulation frameworks (25, 37, 38)"`); 0 `cited.yaml` links. (5 mathematical relations between AFINES and MEDYAN are authored in `popov2016/relations.yaml`.)

---

## 3. `freedman2019 + kierfeld2005`

- **2 verified ModelAtlas relations**:
  1. `rel_kierfeld2005_langmuir_same_as_chem_pot` (`same_as`, symbolic `verified`, `too_plain = False`): `kierfeld2005/linker_site_occupancy` (`X_1 = K*c_x / (1 + K*c_x)`) is symbolically identical to `freedman2019/chemical_potential_from_occupancy` (`N_b = N / (1 + exp(-mu/(k_B T)))`) at `N = 1`, `mu/(k_B T) = ln(K*c_x)`.
  2. `rel_kierfeld2005_lp_2d_scales_as_3d` (`scales_as`, 24 pts, max diff `3.3e-13`): identical `(kappa_B, T)` log-slopes, exposing the exact **2× 2D-vs-3D persistence length discrepancy** — AFINES simulates 2D planar chains (`1` angle per hinge) using the 3D formula `L_p = kappa_B/(k_B T)`, giving an actual 2D tangent-correlation persistence length of `2 * kappa_B / (k_B T) = 34 um` instead of `17 um` (`q_kierfeld2005_2d_vs_3d_persistence_length_and_no_validation`).
- **0 cited links, 0 validations**: neither cites the other, and `kierfeld2005` contains no experimental measurements (`origin: measured`).

---

## 4. `bathe2008 + freedman2019`

- **Honest 0 direct relations, 0 cited links, 0 validations** (`q_bathe2008_no_direct_relation_and_physics_contrast`): Bathe et al. (2008) model continuum macroscopic bending stiffness `kappa_B(N, alpha)` and peristaltic elasticity of homogeneous multi-filament (`N = 4..400`) bundles, disjoint from Freedman et al.'s (2019) 2-filament two-species domain-sorting models. Cross-framework contrasts logged in `questions.yaml` include the reciprocal Langmuir isotherm duality (`delta = delta_min/X` vs `N_b = N*X`), the 8,800× actin extensional-rigidity coarse-graining gap (`E_f A_f = 44 nN` vs AFINES `k_a * l_a = 5 pN`), and Bathe's elastic peristaltic length `lambda_max ~ kappa_f^(1/4)` vs Freedman's rigid-mismatch thermal bending scale `l_g* ~ kappa_f^(1/3)`.

# Audit Log: `weirich2017`

Independent second-agent audit of `corpora/actin_bundling/records/weirich2017/` against `corpus/actin_bundling/2017_Weirich.pdf` (PDF pages 1–16; PNAS 114(9), 2131–2136 + SI Appendix pp. SI-1 to SI-10), inspecting all 16 rendered page images (`pages/page-001.png`–`page-016.png`) and text extractions (`text/page-001.txt`–`page-016.txt`). Figures remain deferred to the digitize stage.

## 1. Fixes Made (Genuine Transcription / Implementation Errors)

- `tactoid_rescaled_free_energy` (p. 8 / SI-2): fixed incorrect small-angle (`alpha < 1e-2`) Taylor series coefficients in `_f_and_g(a)` for `g(alpha) = sin(alpha) - alpha*cos(alpha) - (1/3)*sin^3(alpha)` from `-4/105, +2/405, -16/467775` to the exact expansion `-11/315, +17/3780, -461/1247400` (eliminated `2.3e-6` discontinuity at `alpha = 1e-2`, reducing error to machine precision `~1e-16`).
- `tactoid_equilibrium_minimizer` (p. 8 / SI-2): fixed the same incorrect small-angle Taylor coefficients in `_f_and_g(a)` for `g(alpha)`; improved `rel_tactoid_min_reduces_to_large_w` (`omega -> oo`) asymptotic convergence residual by 60x, from `1.8e-8` to `3.0e-10`.
- `tactoid_shape_dynamics_ode` (pp. 8–9 / SI-2 to SI-3): fixed the same incorrect small-angle Taylor coefficients in `_f_and_g(a)` for `g(alpha)` and the corresponding wrong derivative coefficients in `_fp_and_gp(a)` for `g'(alpha) = alpha*sin(alpha) - sin^2(alpha)*cos(alpha)` from `-4/15, +2/45, -16/42525` to the exact `-11/45, +17/420, -461/113400` (eliminated `3.3e-6` error at `alpha = 1e-2`).
- `c_filamin_main` (p. 1 / 2131): removed erroneous `from_source: ref-19`; `"Adding 0.26 μM of the F-actin cross-linker filamin (19)"` cites Nakamura et al. (2011) as a general review identifying filamin, whereas `0.26 uM` is Weirich et al.'s own experimental concentration (10 mol % of `2.6 uM` actin).

## 2. Additions Made (Missing Printed Equations, Runs, Comparisons, and Open Questions)

- `capillary_dissipation_power` + `run_cap_dissipation_probe` (p. 9 / SI-3): added model for the printed capillary viscous dissipation scaling equation `D ~ eta * (pi * R_0^2 * lambda) * (v / R_0)^2` (`prefactor: unknown`, `built_from: [capillary_axial_velocity]`, with `symbolic()`) and connected it to `cap_axial_vel_out` (`1.508e-20 W`).
- `run_gamma_I_coal_10_printed_eta` + `cmp_gamma_I_coal_10_printed_eta` (p. 4 / 2134 & p. 11 / SI-5): added reproduction run/comparison pairing Table S1 10 mol % coalescence `gamma_I/eta = 0.11 +/- 0.03 um/s` with printed `eta ~ 3 Pa*s` -> `330 nN/m` vs `~300 nN/m` (previously only asserted in prose in `q_gamma_I_filamin_condition_choice`; intentionally retains the validator cross-condition warning noting 10% vs 5% filamin).
- `run_gamma_I_coal_10_chained_eta` + `cmp_gamma_I_coal_10_chained_eta` (p. 4 / 2134 & p. 11 / SI-5): added chained reproduction run/comparison pairing Table S1 10 mol % coalescence `gamma_I/eta` with recomputed Stokes-Einstein `eta = 3.133 Pa*s` -> `344.7 nN/m` vs `~300 nN/m`.
- `run_enrichment_ratio_entropic` + `cmp_enrichment_ratio_entropic` (p. 2 / 2132): added reproduction run/comparison directly pairing the adjacent p. 2 col 1 entropic nematic threshold `c_actin_entropic_nematic_main` (`~250 uM`, Ref. 23) with `c_actin_dilute_main` (`2.6 uM`) against `c_actin_enrichment_factor` (`100` -> `96.15`, `not reproduced` due to `2.6` vs Fig. 1B's `2.5 uM`).
- `run_frap_diff_scaling_fig1c` + `cmp_frap_diff_scaling_fig1c` (p. 2 / 2132 & p. 7 / SI-1): added reproduction run/comparison of `frap_diffusion_scaling` using Fig. 1C's exact fitted `tau_R_frap_fig1c` (`880 s`, yielding `3.09e-3 um^2/s` at `L_b = 1.65 um`, `prefactor: unknown`).
- `run_ode_coal_10_fig2a_omega` + `cmp_ode_Lfinal_coal_10_fig2a_vs_Leq` (p. 4 / 2134, p. 9 / SI-3, p. 11 / SI-5): added ODE run at Fig. 2A's `omega = 1.4` (`omega_fig2a_10_filamin`) with Table S1 10% coalescence `v_cap` and `L_eq`, computing `tau = 14.33 s` (`tau_linearized = 12.78 s`) vs `12.52 s` at Table S1's `omega = 1.87` and verifying `L_final = 6.040 um` (`reproduced`).
- `cmp_filament_len_1p5_vs_short_range`, `cmp_filament_len_0p5_vs_short_range`, `cmp_filament_len_0p25_vs_short_range` (p. 7 / SI-1 & p. 8 / SI-2): added consistency checks verifying computed filament lengths (`0.18 um`, `0.54 um`, `1.08 um`) against SI p. 8's printed experimental range `filament_length_short_range_si` (`0.18 to 2 um`; all `reproduced`).
- `cmp_tactoid_ar_min_vs_scaling_contr10` (p. 8 / SI-2 & p. 11 / SI-5): added approximation comparison between exact circular-arc minimization (`7.751`) and `2*omega^(1/2)` (`7.746`) at Table S1 10% contraction `omega = 15.0`.
- `cmp_strain_contr_5_vs_range` (p. 5 / 2135 & p. 11 / SI-5): added consistency check for 5 mol % contraction strain (`34.8%` vs `20 to 60%`, `reproduced`).
- `cmp_ode_Lfinal_coal_5_vs_Leq`, `cmp_ode_Lfinal_contr_2p5_vs_Leq`, `cmp_ode_Lfinal_contr_5_vs_Leq` (p. 9 / SI-3 & p. 11 / SI-5): added equilibrium-relaxation consistency checks across the remaining three Table S1 ODE runs (all `reproduced`).
- `q_fig3d_omega_vs_volume_separation` (p. 4 / 2134 & p. 9 / SI-3): logged new `unexplained_difference` finding that at fixed volume `V`, changing `omega` from `2.0` to `1.4` shifts `tau` by only `~5%`, whereas Fig. 3D's dashed curves are separated by `~2.7-3x` (tracking Table S1's `2.2x` difference in `L_eq`, `2.74` vs `6.04 um`, rather than `omega` alone).
- `q_inverted_eta_over_gamma_text_vs_fig3d_tableS1` (p. 4 / 2134, p. 9 / SI-3, p. 11 / SI-5): logged new `possible_error` finding that Main Text p. 4 col 1 and SI p. 9 col 1 write `eta/gamma_I` when describing Fig. 3D and Table S1, while both Fig. 3D (`gamma_I/eta (nm/s)`) and Table S1 (`gamma_I/eta, um/s`) plot and tabulate the reciprocal `gamma_I/eta`.
- `q_cross_cp_condition_mixing_and_movies3_mismatch` (p. 2 / 2132, p. 3 / 2133, p. 4 / 2134, p. 11 / SI-5, p. 12 / SI-6): logged new `unexplained_difference` finding documenting (1) Fig. 3A (`1.5 mol % CP`) vs Movie S3 (`1 mol % CP`) mismatch for coarsening, and (2) p. 4 combining `1 mol % CP` (`270 nm`) coalescence fits (Figs. 3B–C, Table S1) with `1.5 mol % CP` (`180 nm`) Fig. 2A `omega` and Fig. 1C FRAP viscosity without addressing the 50% filament-length difference.
- `q_isotropic_dEdt_sign_and_capillary_velocity_pi` (p. 8 / SI-2 & p. 9 / SI-3): logged new `possible_error` finding documenting (1) SI p. 8 `dE/dt ~ -gamma_I * r * dL/dt` sign mismatch for shortening (`dL/dt < 0`), and (2) SI p. 9 equality `v = 2*pi*(du/dt)/k` vs exact linearized incompressible continuity amplitude `v_z = (2/k)*(du/dt)` (factor of `pi` from order-of-magnitude `d(v_z)/dz ~ v/lambda` scaling).

## 3. Disagreements Logged as `kind: reading` Questions (First-Reader Judgements Preserved)

- `q_audit_reading_tableS1_coal_conditions` (p. 11 / SI-5): disagreed with asymmetric condition annotation in `data.yaml`, where `capping_protein_mole_fraction: "0.25 percent"` is set on `tableS1_contr_*` (from Figs. 4B–D) while `capping_protein_mole_fraction: "1 percent"` is omitted on `tableS1_coal_5_*` and `tableS1_coal_10_*` (from Figs. 3B–C and Movie S4), hiding the 1% vs 1.5% CP mismatch in cross-figure runs; preserved first-reader choice and logged as `reading`.
- `q_audit_reading_ambient_temp_and_frap_spot` (p. 7 / SI-1 & p. 2 / 2132): disagreed with (a) using vesicle-preparation temperature `ambient_temperature_si` (`~20 degC`, p. 7) directly as the FRAP imaging temperature instead of an `added.yaml` input (minimal numerical impact: `<2%` in Kelvin, `eta = 3.13` vs `3.19 Pa*s` at 25 °C), and (b) setting `frap_bleached_size_added` baseline to `"1.65"` um (equal to `sqrt(D*tau_R) = 1.643 um`, mildly back-calculated despite `prefactor: unknown`, though alternatives `["1.0", "2.0", "2.5"]` span Fig. 1C); preserved first-reader choices and logged as `reading`.
- `q_audit_reading_mesh_and_ratio_from_source` (p. 1 / 2131 & p. 4 / 2134): disagreed with `from_source: ref-18` on `mesh_size_dilute_main` (Ref. 18 supplies the formula `xi = 0.3 c^(-1/2)`, not the numerical value `~1 um` for Weirich et al.'s `2.6 uM` sample) and `from_source: ref-24` on `gamma_I_ratio_protein_droplets` (`~10` is Weirich et al.'s own comparison ratio against both Refs. 24 and 29); preserved first-reader choices and logged as `reading`.

## 4. Verification of `data.yaml` Against Page Images (All 73 Entries)

- `c_actin_dilute_main` (p. 1 / 2131): verified `"2.6"` `uM`, origin `reported`, quote matches `page-001.png` (`2.6 μM`).
- `c_actin_dilute_fig1b` (p. 2 / 2132): verified `"≈ 2.5"` `uM`, origin `reported`, matches Fig. 1B on `page-002.png`.
- `c_actin_monomeric_si` (p. 7 / SI-1): verified `"2.64"` `uM`, origin `reported`, matches `page-007.png`.
- `c_actin_tmr_si` (p. 7 / SI-1): verified `"0.32"` `uM`, origin `reported`, matches `page-007.png`.
- `c_actin_tactoid_est` (p. 1 / 2131): verified `"250"` `uM`, origin `estimated`, matches `page-001.png` (`250-μM`).
- `c_actin_entropic_nematic_main` (p. 2 / 2132): verified `"~250"` `uM`, `from_source: ref-23`, matches `page-002.png`.
- `c_actin_entropic_nematic_fig1b` (p. 2 / 2132): verified `"≥ 250"` `uM`, matches Fig. 1B on `page-002.png`.
- `c_actin_enrichment_factor` (p. 2 / 2132): verified `"100"` `dimensionless`, origin `derived`, matches `page-002.png` (`100-fold`).
- `c_filamin_main` (p. 1 / 2131): verified `"0.26"` `uM` on `page-001.png`; fixed erroneous `from_source: ref-19` (see Section 1).
- `filamin_frac_fig1a` (p. 2 / 2132): verified `"10"` `percent`, matches Fig. 1A caption on `page-002.png`.
- `filamin_frac_tactoid_range` (p. 2 / 2132): verified `"1 to 10"` `percent`, matches `page-002.png`.
- `filamin_frac_fig2a_range` (p. 3 / 2133): verified `"2.5 to 15"` `percent`, matches `page-003.png`.
- `filamin_frac_2p5` (p. 3 / 2133): verified `"2.5"` `percent`, matches `page-003.png`.
- `filamin_frac_5` (p. 2 / 2132): verified `"5"` `percent`, matches Fig. 1C caption on `page-002.png`.
- `filamin_frac_10` (p. 4 / 2134): verified `"10"` `percent`, matches Fig. 3B caption on `page-004.png`.
- `filamin_frac_15` (p. 3 / 2133): verified `"15"` `percent`, matches `page-003.png`.
- `cp_frac_1p5` (p. 2 / 2132): verified `"1.5"` `percent`, matches Fig. 1C caption on `page-002.png`.
- `cp_frac_1p0` (p. 2 / 2132): verified `"1"` `percent`, matches Fig. 1A caption on `page-002.png`.
- `cp_frac_0p5` (p. 5 / 2135): verified `"0.5"` `percent`, matches Fig. 4A caption on `page-005.png`.
- `cp_frac_0p25` (p. 5 / 2135): verified `"0.25"` `percent`, matches Fig. 4B caption on `page-005.png`.
- `actin_monomer_spacing_si` (p. 7 / SI-1): verified `"~2.7"` `nm`, matches `page-007.png` (`∼1 monomer per 2.7 nm`).
- `filament_length_avg_main` (p. 1 / 2131): verified `"~180"` `nm`, matches `page-001.png`.
- `filament_length_tactoid_range` (p. 2 / 2132): verified `"0.1 to 1"` `um`, matches `page-002.png`.
- `filament_length_short_range_si` (p. 8 / SI-2): verified `"0.18 to 2"` `um` (with `note` documenting mixed-unit print `between 180 nm and 2 μm` on `page-008.png`).
- `filament_length_homogeneous_si` (p. 8 / SI-2): verified `"~200"` `nm`, matches `page-008.png`.
- `filament_length_bipolar_ref23` (p. 8 / SI-2): verified `"~200"` `nm`, `from_source: ref-23`, matches `page-008.png`.
- `filament_persistence_length_si` (p. 8 / SI-2): verified `"~10"` `um`, matches `page-008.png`.
- `filamin_length_main` (p. 3 / 2133): verified `"~150"` `nm`, `from_source: ref-19`, matches `page-003.png`.
- `mesh_size_dilute_main` (p. 1 / 2131): verified `"~1"` `um` on `page-001.png`; logged `from_source: ref-18` reading disagreement (`q_audit_reading_mesh_and_ratio_from_source`).
- `confocal_depth_si` (p. 7 / SI-1): verified `"~1"` `um`, matches `page-007.png`.
- `ambient_temperature_si` (p. 7 / SI-1): verified `"~20"` `degC` on `page-007.png`; logged vesicle-prep vs imaging temperature reading disagreement (`q_audit_reading_ambient_temp_and_frap_spot`).
- `tau_R_frap_fig1c` (p. 2 / 2132): verified `"880"` `s`, origin `fitted`, matches Fig. 1C caption on `page-002.png`.
- `tau_R_frap_text` (p. 2 / 2132): verified `"~900"` `s`, origin `fitted`, matches `page-002.png` (and p. 7 / SI-1).
- `D_frap_est` (p. 2 / 2132): verified `"~0.3 x 10^-2"` `um**2/s`, matches `page-002.png` (and p. 7 / SI-1).
- `filament_half_length_frap_si` (p. 8 / SI-2): verified `"~90"` `nm`, matches `page-008.png` (`length L ∼ 90 nm` with p. 7 axes `2L` and `2r`).
- `filament_diameter_frap_si` (p. 8 / SI-2): verified `"~7"` `nm`, matches `page-008.png` (`diameter 2r ∼ 7 nm`).
- `eta_frap_est` (p. 8 / SI-2): verified `"~3"` `Pa*s`, matches `page-008.png` (and p. 2 / 2132).
- `gamma_I_est_main` (p. 4 / 2134): verified `"~300"` `nN/m`, matches `page-004.png`.
- `gamma_I_ratio_protein_droplets` (p. 4 / 2134): verified `"~10"` `dimensionless` on `page-004.png`; logged single-ref `from_source: ref-24` reading disagreement (`q_audit_reading_mesh_and_ratio_from_source`).
- `tactoid_length_growth_60min` (p. 1 / 2131): verified `"~1 to 4"` `um`, matches `page-001.png`.
- `tactoid_length_homogeneous_si` (p. 8 / SI-2): verified `"~5"` `um`, matches `page-008.png`.
- `tactoid_length_bipolar_ref23` (p. 8 / SI-2): verified `"~50"` `um`, `from_source: ref-23`, matches `page-008.png`.
- `aspect_ratio_2p5_filamin` (p. 3 / 2133): verified `"~3"` `dimensionless`, matches `page-003.png`.
- `aspect_ratio_15_filamin` (p. 3 / 2133): verified `"~2"` `dimensionless`, matches `page-003.png`.
- `omega_range_fig2a_si` (p. 8 / SI-2): verified `"1.1 to 2.2"` `dimensionless`, matches `page-008.png` (`1.1 ≤ ω ≤ 2.2`).
- `omega_fig2a_5_filamin` (p. 4 / 2134): verified `"2"` `dimensionless`, matches `page-004.png` (`ω = 2`; Fig. 3D prints `ω=2.0`).
- `omega_fig2a_10_filamin` (p. 4 / 2134): verified `"1.4"` `dimensionless`, matches `page-004.png`.
- `coarsening_exponent_measured` (p. 3 / 2133): verified `"0.47 ± 0.01"` `dimensionless`, matches `page-003.png` (and Fig. 3A caption p. 4).
- `coarsening_exponent_theory` (p. 3 / 2133): verified `"0.5"` `dimensionless`, `from_source: ref-28`, matches `page-003.png`.
- `tableS1_coal_5_omega` (p. 11 / SI-5): verified `"1.90 ± 0.53"` `dimensionless` on `page-011.png`; logged CP-condition asymmetry (`q_audit_reading_tableS1_coal_conditions`).
- `tableS1_coal_5_gamma_eta` (p. 11 / SI-5): verified `"0.14 ± 0.02"` `um/s` on `page-011.png`.
- `tableS1_coal_5_Leq` (p. 11 / SI-5): verified `"2.74 ± 0.52"` `um` on `page-011.png`.
- `tableS1_coal_10_omega` (p. 11 / SI-5): verified `"1.87 ± 0.07"` `dimensionless` on `page-011.png`.
- `tableS1_coal_10_gamma_eta` (p. 11 / SI-5): verified `"0.11 ± 0.03"` `um/s` on `page-011.png`.
- `tableS1_coal_10_Leq` (p. 11 / SI-5): verified `"6.04 ± 0.38"` `um` on `page-011.png`.
- `tableS1_contr_2p5_omega` (p. 11 / SI-5): verified `"23.9 ± 14.4"` `dimensionless` on `page-011.png`.
- `tableS1_contr_2p5_gamma_eta` (p. 11 / SI-5): verified `"2.5 ± 0.9 x 10^-3"` `um/s` on `page-011.png` (and confirmed `boxgraph.precision.parse` yields `0.0025 +/- 0.0009`).
- `tableS1_contr_2p5_Leq` (p. 11 / SI-5): verified `"5.77 ± 0.55"` `um` on `page-011.png`.
- `tableS1_contr_5_omega` (p. 11 / SI-5): verified `"11.9 ± 5.48"` `dimensionless` on `page-011.png`.
- `tableS1_contr_5_gamma_eta` (p. 11 / SI-5): verified `"4.5 ± 0.7 x 10^-3"` `um/s` on `page-011.png` (`0.0045 +/- 0.0007`).
- `tableS1_contr_5_Leq` (p. 11 / SI-5): verified `"5.87 ± 2.01"` `um` on `page-011.png`.
- `tableS1_contr_10_omega` (p. 11 / SI-5): verified `"15.0 ± 7.78"` `dimensionless` on `page-011.png`.
- `tableS1_contr_10_gamma_eta` (p. 11 / SI-5): verified `"3.5 ± 0.7 x 10^-3"` `um/s` on `page-011.png` (`0.0035 +/- 0.0007`).
- `tableS1_contr_10_Leq` (p. 11 / SI-5): verified `"5.68 ± 0.54"` `um` on `page-011.png`.
- `bundle_initial_length` (p. 4 / 2134): verified `"> 50"` `um`, matches `page-004.png`.
- `bundle_initial_radius_main` (p. 5 / 2135): verified `"~200"` `nm`, matches `page-005.png`.
- `bundle_initial_radius_si` (p. 9 / SI-3): verified `"~200"` `nm`, matches `page-009.png`.
- `bundle_omega_capillary_si` (p. 9 / SI-3): verified `"~2"` `dimensionless`, matches `page-009.png`.
- `capillary_lambda_max_predicted` (p. 9 / SI-3): verified `"~4"` `um`, matches `page-009.png`.
- `capillary_lambda_observed` (p. 9 / SI-3): verified `"~10"` `um`, matches `page-009.png`.
- `capillary_lambda_to_R0_observed` (p. 9 / SI-3): verified `"~10 to 100"` `dimensionless`, matches `page-009.png`.
- `bundle_shape_change_timescale` (p. 5 / 2135): verified `"~10"` `min`, matches `page-005.png`.
- `bundle_contractile_strain_range` (p. 5 / 2135): verified `"20 to 60"` `percent`, matches `page-005.png`.

## 5. Verification of Models Against Printed Equations, Units, and `symbolic()` (All 46 Models)

Every closed-form model was verified against the page images (`page-001.png`–`page-010.png`), `pint` dimensional conversions, and `run()` vs `symbolic()` agreement (`<= 2.0e-15` relative difference).

- `actin_filament_length_cp` (p. 7 / SI-1): verified `N_m = 1/r_CP`, `L_f = a_m/r_CP`, `from_words: true`, matches text on `page-007.png`.
- `actin_mass_concentration` (p. 7 / SI-1): verified `c [mg/mL] = c_A [uM] * 1e-6 * M_w [g/mol]` via `pint`.
- `actin_mesh_size` (p. 7 / SI-1): verified `xi = 0.3 * c^(-1/2)` (`described_from: ref-18`) against `page-007.png`.
- `anchoring_strength_definition` (p. 3 / 2133): verified `omega = gamma_A / gamma_I` against `page-003.png`.
- `anchoring_strength_from_aspect_ratio` (p. 3 / 2133): verified `omega = (L/r)^2 / 4` (inversion of `L/r = 2*omega^(1/2)`) against `page-003.png`.
- `bundle_contractile_strain` (p. 5 / 2135): verified `strain = (L_i - L_f) / L_i` against `page-005.png` & Fig. 4E.
- `capillary_axial_velocity` (p. 9 / SI-3): verified printed equality `v = 2*pi*(du/dt)/k` against `page-009.png`; logged continuity factor-of-`pi` note in `q_isotropic_dEdt_sign_and_capillary_velocity_pi`.
- `capillary_dissipation_power` (p. 9 / SI-3): added during audit for printed `D ~ eta*(pi*R_0^2*lambda)*(v/R_0)^2`; verified units (`1e-18` W) and `run()` vs `symbolic()`.
- `capillary_dispersion_printed` (p. 9 / SI-3): verified printed `Omega(k) ~ (gamma_I/(eta*R_0))*k^2*[1 - (1+omega)*k^2]` (`prefactor: unknown`) against `page-009.png`.
- `capillary_dispersion_corrected` (p. 9 / SI-3): verified evidenced `variant_of` replacing `(1+omega)` with `(1+2*omega)` from SymPy O(`u^2`) expansion of printed integrals.
- `capillary_dynamic_factor_scaling` (p. 9 / SI-3): verified `Lambda(k) ~ (gamma_I/(eta*R_0))*k^2` (`prefactor: unknown`) against `page-009.png`.
- `capillary_energy_cost_printed` (p. 9 / SI-3): verified printed `delta_F = (pi*lambda*u^2*gamma_I/(2*R_0))*[(1+omega)*k^2 - 1]` and unit factor `1e-21 J/(nN/m * um^2)` against `page-009.png`.
- `capillary_energy_cost_corrected` (p. 9 / SI-3): verified evidenced `variant_of` with `(1+2*omega)*k^2 - 1` via SymPy expansion of printed `A`, `A_omega`, `F`.
- `capillary_fastest_mode_printed` (p. 9 / SI-3): verified printed `k_max = 1/sqrt(2*(1+omega))` (`lambda_max = 3.078 um` at `R_0 = 200 nm, omega = 2`) against `page-009.png`.
- `capillary_fastest_mode_corrected` (p. 9 / SI-3): verified evidenced `variant_of` `k_max = 1/sqrt(2*(1+2*omega))` (`lambda_max = 3.974 um` matching printed `~4 um`).
- `capillary_free_energy_integral` (p. 9 / SI-3): verified exact quadrature of printed `A`, `A_omega`, `F`, cancellation-free `dA`, and periodic trapezoidal convergence (`ref_capillary_integral_grid`: `2.8e-14`).
- `capillary_instability_threshold_printed` (p. 9 / SI-3): verified printed `lambda_c = 2*pi*R_0*sqrt(1+omega)` against `page-009.png` & Fig. S3.
- `capillary_instability_threshold_corrected` (p. 9 / SI-3): verified evidenced `variant_of` `lambda_c = 2*pi*R_0*sqrt(1+2*omega)`.
- `capillary_rate_scale` (p. 9 / SI-3): verified `rate_0 = gamma_I / (eta * R_0)` (`nN/m / (Pa*s * nm) = 1/s`) via `pint`.
- `capillary_volume_conservation_printed` (p. 9 / SI-3): verified printed `2*pi*int_0^lambda rho^2 dz = 2*pi*lambda*(R_0^2 + u^4/(16*R_0^2))` against `page-009.png` via SymPy.
- `capillary_volume_conservation_corrected` (p. 9 / SI-3): verified evidenced `variant_of` with `pi*int_0^lambda rho^2 dz`.
- `capillary_wavelength_from_wavenumber` (p. 10 / SI-4): verified `lambda = 2*pi*R_0/k` against Fig. S3 on `page-010.png`.
- `capillary_wavenumber_definition` (p. 9 / SI-3): verified `k = 2*pi*R_0/lambda` against `page-009.png`.
- `coalescence_tau_linear_scaling` (p. 4 / 2134): verified `slope = eta/gamma_eff`, `tau = slope*(L_i - L_f)` and `1e3` `s/um` unit factor against `page-004.png`.
- `concentration_enrichment_ratio` (p. 2 / 2132): verified `ratio = c_dense / c_dilute` against `page-002.png`.
- `director_crossover_parameter` (p. 8 / SI-2): verified `chi = gamma_A * V^(1/3) / K` (`1e-3` factor for `nN/m * um / pN`) against `page-008.png`.
- `filament_half_length_from_length` (p. 7 / SI-1): verified `L = L_f / 2` from axes `2L` and `2r` on `page-007.png`.
- `filament_radius_from_diameter` (p. 8 / SI-2): verified `r = d / 2` from `2r ~ 7 nm` on `page-008.png`.
- `filamin_concentration_from_fraction` (p. 1 / 2131): verified `c_fil = c_A * f_fil` against `page-001.png` & Fig. 1A.
- `frap_diffusion_scaling` (p. 7 / SI-1): verified `D ~ L_b^2 / tau_R` (`from_words: true`, `prefactor: unknown`) against `page-007.png`.
- `interfacial_tension_from_capillary_velocity` (p. 4 / 2134): verified `gamma_I [nN/m] = 1e3 * eta [Pa*s] * v_cap [um/s]` via `pint`.
- `isotropic_coalescence_tau_energy_balance` (p. 8 / SI-2): verified `tau ~ (eta/gamma_I)*(V/(r*L_i^2))*(L_i - L_f)` (`prefactor: unknown`) against `page-008.png`.
- `isotropic_relaxation_rate_printed` (p. 8 / SI-2): verified faithful implementation of dimensionally inverted printed `dL/dt ~ -L*(r*eta/gamma_I)`, `tau ~ gamma_I/eta` (`prefactor: unknown`) against `page-008.png`.
- `isotropic_relaxation_rate_corrected` (p. 8 / SI-2): verified evidenced `variant_of` `dL/dt ~ -L*(gamma_I/(r*eta))`, `tau ~ r*eta/gamma_I` from balancing printed `dE/dt` and `D`.
- `stokes_einstein_ellipsoid_diffusion` (p. 8 / SI-2): verified printed `D = k_B*T*ln(2*L/r)/(6*pi*eta*L)` (`described_from: ref-43`) and `pint` conversion against `page-008.png`.
- `stokes_einstein_ellipsoid_viscosity` (p. 8 / SI-2): verified exact inversion `eta = k_B*T*ln(2*L/r)/(6*pi*D*L)` (`3.133 Pa*s` at `20 degC, L=90 nm, r=3.5 nm, D=0.003 um^2/s`).
- `tactoid_aspect_ratio_large_omega` (p. 3 / 2133): verified `L/r = 2*omega^(1/2)` (`described_from: ref-20`) against `page-003.png`.
- `tactoid_aspect_ratio_scaling` (p. 8 / SI-2): verified piecewise `2*omega^(1/2)` (`omega >= 1`), `1 + omega` (`0 <= omega < 1`) against `page-008.png`.
- `tactoid_coarsening_power_law` (p. 3 / 2133): verified `L ~ t^alpha` (`role: fit`, `prefactor: unknown`) against `page-003.png`.
- `tactoid_dissipation_power` (p. 9 / SI-3): verified printed `D = (1/2)*V*[eta_R*((1/R)*dR/dt)^2 + eta_alpha*((1/alpha)*dalpha/dt)^2]` and `1e-18 W` unit factor against `page-009.png`.
- `tactoid_equilibrium_minimizer` (p. 8 / SI-2): verified minimization of `F_hat(alpha, omega)` and `L/r = cot(alpha_eq/2)`, plus fixed `g(alpha)` small-angle Taylor series (Section 1); independently confirmed circular-arc asymptotics `L/r -> 2*omega^(1/2)` (`omega -> oo`) and `1 + 2*pi/(3*(pi^2 - 8))*omega = 1 + 1.1202*omega` (`omega -> 0`).
- `tactoid_exponential_relaxation` (p. 3 / 2133): verified printed `L(t) = L_f + (L_i - L_f)*exp(-t/tau)` against `page-003.png`.
- `tactoid_free_energy` (p. 8 / SI-2): verified printed `F = gamma_I*A + gamma_A*A_omega` and `1e-21 J` conversion against `page-008.png`.
- `tactoid_geometry_circular_arc` (p. 8 / SI-2 & p. 10 / SI-4): independently derived all circular-arc integrals via SymPy and confirmed exact identities `A = 4*pi*R^2*f`, `A_omega = 2*pi*R^2*g`, `V = 2*pi*R^3*g = R*A_omega`, `L = 2*R*sin(alpha)`, `r = 2*R*(1 - cos(alpha))`, matching `page-008.png` and Fig. S2B on `page-010.png`.
- `tactoid_rescaled_free_energy` (p. 8 / SI-2 & p. 10 / SI-4): verified `F/(gamma_I*V^(2/3))`, fixed small-angle `g(alpha)` Taylor series (Section 1), and independently confirmed Fig. S2C (`page-010.png`) is rescaled downward by `(2*pi)^(1/3) = 1.8453` (`q_figS2C_missing_2pi_one_third_factor`).
- `tactoid_shape_dynamics_ode` (p. 9 / SI-3): independently derived `dF/dt + D = 0` with `dV/dt = 0` via SymPy confirming exact `dF_hat/dalpha` and mobility factor, verified RK4 step-size convergence across all 5 Table S1 runs (`rel_diff <= 8.7e-7` from `n_steps = 400` to `3200`), and fixed small-angle Taylor series for `g` and `g'` (Section 1).

## 6. Verification of `added.yaml`, `relations.yaml`, `refinements.yaml`, `references.yaml`, and `source.yaml`

- `added.yaml` (all 23 entries, pp. 1–11): verified every entry has a physical or figure-consistent basis and real non-trivial alternatives (no exact mathematical identities given spurious alternatives); noted `frap_bleached_size_added` baseline choice in `q_audit_reading_ambient_temp_and_frap_spot` and confirmed `contr_2p5_Li_added` alternative `"7.2"` genuinely flips `cmp_strain_contr_2p5_vs_range` from `reproduced` to `close` (`19.86%`).
- `relations.yaml` (all 11 entries, pp. 3–9): ran `boxgraph relations --only weirich2017`; 10/11 verified (including `rel_tactoid_min_reduces_to_large_w` improving from `1.8e-8` to `3.0e-10` after the Taylor fix and `rel_cap_integral_reduces_to_energy_corrected` verified to `1.7e-8`) while `rel_cap_integral_reduces_to_energy_printed` is honestly `not_verified` (`0.9387` residual), exposing the paper's missing factor of 2 on `omega`.
- `refinements.yaml` (all 6 entries, pp. 8–9): verified all 6 numerical convergence series settle cleanly (`<= 2.0e-8` shifts on `xtol`, `n_bracket`, `n_steps` for `tau` and `alpha_min`, and `n_grid`).
- `references.yaml` (all 49 entries, p. 6 / 2136): verified labels `1`–`49`, author names, years, and cross-references match `page-006.png`.
- `source.yaml` (pp. 1–16): verified title, authors, year, venue, DOI, page mapping, and notes against `page-001.png`–`page-016.png` (confirming Fig. 4E y-axis tick misprint `0.2, 0.6, 0.4, 0.5, 0.6` on `page-005.png`, SI p. 7 `~2 uM` typo on `page-007.png`, and SI figure/movie cross-reference typos on pp. 7–9).

## 7. Final Validation Result

- Command: `.venv/bin/boxgraph validate corpora/actin_bundling/records/weirich2017 --cache build/cache`
- Result: `weirich2017: SOUND (0 errors, 2 warnings)` — the 2 warnings are intentional cross-filamin-concentration condition notices on `run_gamma_I_coal_10_printed_eta` and `run_gamma_I_coal_10_chained_eta` (`10 percent | 5 percent`), explicitly documenting the paper's pairing of its single 5%-filamin FRAP viscosity with Table S1's 10%-filamin coalescence velocity.

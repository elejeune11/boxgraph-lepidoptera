# Audit Log — `chandrasekaran2024`

Independent second-agent audit of `corpora/actin_bundling/records/chandrasekaran2024/` against
`corpus/actin_bundling/2024_Chandrasekaran.pdf` (PDF pages 1–37, including main text, Methods,
Supplementary Information pp. 15–36, and Supplementary Movie descriptions p. 37), verified directly
against rendered page images (`pages/page-NNN.png`) and extracted text (`text/page-NNN.txt`).

- **Final validation status:** `chandrasekaran2024: SOUND (0 errors, 2 warnings)` (both warnings are
  for the two `runnable: false` stochastic many-body models `cytosim_actin_vasp_droplet` and
  `lammps_vasp_actin_droplet`).
- **Scope counts after audit:** 118 data entries, 46 models (39 runnable printed/from-words + 2
  corrected variants + 1 reconstructed Stokes-Einstein check + 4 non-runnable many-body/external
  models; 42 runnable total), 14 added inputs, 47 runs, 23 comparisons, 3 numerical refinements,
  4 verified stated relations, 27 open questions.

---

## 1. Fixes Made (Genuine Transcription Errors)

- `vasp_polymerization_speedup` (p. 5): Fixed `printed` from `"3"` to `"≤ 3"` to preserve the
  printed upper-bound qualifier (`"up to three times faster than free barbed ends"`, `page-005.png`),
  consistent with `exp_gabmut_disk_pct_1_7p5` (`"< 10"`) and `vasp_disordered_pct` (`"about 60"`).

---

## 2. Additions Made (Unrecorded Printed Parameters, Equations, Runs, Comparisons, Questions)

### Added printed data (`data.yaml`)
- `lammps_filament_count_low` (p. 21): Added printed LAMMPS low filament count `"6"` (`page-021.png`).
- `lammps_filament_count_high` (p. 21): Added printed LAMMPS high filament count `"20"` (`page-021.png`).
- `lammps_filaments_added_increment` (p. 21): Added printed LAMMPS equilibration filament addition increment `"2"` (`page-021.png`).
- `suppfig11_filament_counts_range` (pp. 31–32): Added printed Supplementary Fig. 11A, D filament count sweep `"30 to 165"` (`{30, 75, 120, 165}`, `page-031.png` & `page-032.png`).
- `supptable3_ring_pct_1_0p1` (p. 19): Added second probabilistic cell `(k_bind=1.0 1/s, k_unbind=0.1 1/s)` Ring abundance `"66.67"`% from Supplementary Table 3 image (`page-019.png`).
- `supptable3_shell_sb_pct_1_0p1` (p. 19): Added second probabilistic cell `(k_bind=1.0 1/s, k_unbind=0.1 1/s)` Shell with strong bundling abundance `"33.33"`% from Supplementary Table 3 image (`page-019.png`).

### Added models, added inputs, runs, comparisons, and relations
- `lammps_reduced_temp_from_lj_eps` (pp. 16, 21): Added closed-form model (with `symbolic()`) evaluating printed `T* = k_B T / eps_{V-V}` directly in the `k_B T` LJ energy units of Supplementary Table 1 (`T_star = 1 / eps_VV_lj`).
- `run_lammps_reduced_temp_vs_table1` (pp. 16, 21): Evaluated `lammps_reduced_temp_from_lj_eps` on printed `lammps_eps_vasp_vasp = 1.66`, yielding `0.6024` with zero added inputs.
- `cmp_lammps_reduced_temp_vs_table1` (pp. 16, 21): Added non-circular `consistency` comparison of `0.6024` against printed `lammps_reduced_temp = 1.0` -> verdict `not reproduced`.
- `rel_lammps_reduced_temp_lj` (p. 21): Added stated `same_as` relation linking `lammps_reduced_temperature` to `lammps_reduced_temp_from_lj_eps` under `eps_VV_lj = eps_VV / kBT` (verified by `boxgraph relations`).
- `lammps_radial_number_density` (p. 22): Added processing model (with `symbolic()`) for the explicitly printed y-axis equation `Radial density (N/(4*pi*r^2*dr))` in Supplementary Fig. 1B (`page-022.png`).
- `lammps_shell_count_probe_added` (p. 22): Added probe shell count `3000` (alternatives `2500`, `3500`) to run `lammps_radial_number_density` (`run_lammps_radial_density_probe` -> `0.780 / sigma^3`, matching the droplet interior plateau in Supplementary Fig. 1B).

### Added open questions (`questions.yaml`)
- `q_lammps_tstar_vs_eps_vv_contradiction` (pp. 16, 21): `kind: possible_error` — contradiction between `T* = k_B T / eps_{V-V} = 1.0` (p. 21) and `eps_{V-V} = 1.66` in `k_B T` LJ units (p. 16, citing Cho & Kim 2012 Fig. 3 where phase separation requires `k_B T / eps_{V-V} = 0.602 < T_c* ~ 1.08`).
- `q_fig4_text_vs_panel_length_regimes` (p. 5): `kind: possible_error` — main text maps `L_fil^max < R_drop`, `R_drop < L_fil^max < 2 R_drop`, and `L_fil^max > 2 R_drop` to Fig. 4A `(i-ii)`, `(iii-iv)`, `(v-viii)`, contradicting the actual lengths (`ii = pi*R/3 = 1.047 R > R`, `iv = 2*pi*R/3 = 2.094 R > 2 R`) and the regime brackets printed directly beneath Fig. 4A (`i`, `ii-iii`, `iv-viii`).
- `q_fig2_fig3_axis_column_and_supptable3_1_10_errors` (pp. 3, 4, 19): `kind: possible_error` — (1) Supplementary Table 3 classifies `(k_bind = 1.0, k_unbind = 10)` as `'Shell with strong bundling'`, contradicting p. 3 text, Fig. 2, Fig. 3B, and Fig. 3D (>95% unbound VASP, zero bound to >=2 filaments); (2) p. 3 text & Fig. 2 caption swap Fig. 2 rows (`k_bind`) and columns (`k_unbind`); (3) p. 4 text refers to `'Fig. 3D, left column'` and `'Fig. 3D right column'` despite Fig. 3D being a single vertical stack with unbound shells on the right (`k_unbind = 10 s^-1`).
- `q_lammps_lj_shift_and_actin_cutoff_ambiguity` (pp. 16, 20): `kind: possible_error` — (1) unshifted LJ formula printed on p. 20 (with `"Crowder molecules interact through an LJ potential U_{V-V}(r)"` typo) omits whether WCA/LJ cutoffs are shifted to zero; (2) `Actin-Actin` cutoff `1.5 sigma = 15 * sigma_{A-A} >> 2^(1/6) * 0.1 sigma` retains a full `-1.0 k_B T` attractive LJ well between actin beads.
- `q_minor_crossref_and_legend_typos` (pp. 4, 6, 10, 11, 18, 19, 31): `kind: possible_error` — 7 additional verifiable cross-reference/label typos (Fig. 3E cites Fig. 1 for Fig. 2; Fig. 5C caption writes `'rings with strong bundling'`; p. 6 omits `'3E'` after `'Figure'`; Methods p. 10 cites Fig. 1A/1B instead of Fig. 1B/1D; Methods p. 11 prints `'often replicates'` for `'of ten replicates'`; Supp. Table 2 pp. 18–19 cites Figure 4 instead of Figure 2; Supp. Fig. 11D legend writes `'Ring-Shell w. strong bundling'` vs Panel C `'Shell w. strong bundling'`).

---

## 3. Disagreements Logged (First-Reader Judgements Preserved)

- `q_audit_lammps_eps_si_circular_reading` (pp. 16, 21): `kind: reading` — Disagreement on `lammps_eps_si_added` (`4.11 pN*nm`) and `cmp_lammps_reduced_temp`. Back-solving `eps_{V-V} = 4.11 pN*nm` from `T* = 1.0` makes `cmp_lammps_reduced_temp` (`4.11 / 4.11 = 1.0`) circular, whereas the Supplementary Table 1 reading `eps_{V-V} = 1.66 k_B T = 6.82 pN*nm` (listed as an alternative) flips it to `not reproduced` (`0.6026`), matching our non-circular `cmp_lammps_reduced_temp_vs_table1`. First-reader `lammps_eps_si_added` preserved unchanged.
- `q_audit_pca_communality_standardization_reading` (pp. 25, 30, 31, 33): `kind: reading` — Nuance on first-reader note in `pca_var_suppfig4_caption` / `q_pca_variance_pc_count_contradictions` equating mean squared Varimax loading (`sum L^2 / 5`) with PCA explained variance. Checking all four printed Varimax matrices gives `95.74%` (Supp. Fig. 4B, 3 PCs vs `95.83%` caption), `94.43%` (Supp. Fig. 12B, 3 PCs vs `94.72%` caption), `98.67%` (Supp. Fig. 11B.ii, 3 PCs vs `96.86%` caption), and `85.64%` (Supp. Fig. 9B, 2 PCs vs `~95%` 2-PC bar and `99.06%` caption). The `85.64%` vs `~95%` gap in Supp. Fig. 9B shows unweighted mean squared loading only equals PCA explained variance for standardized (correlation-matrix) PCA, not unstandardized covariance-matrix PCA, while still confirming Supp. Figs. 4B and 11B.ii retain 3 PCs.

---

## 4. Verification of Metadata, References, and `data.yaml` vs Page Images

### Metadata & References
- `source.yaml` (pp. 1–37): Checked title, authors, year (2024), venue (`Nature Communications 15, 3139`), DOI (`10.1038/s41467-024-46726-6`), and page label map (`1`–`14`, `SI-1`–`SI-22`, `Movies`) against `page-001.png`–`page-037.png` — exact match.
- `references.yaml` (pp. 11–13, 35–36): Verified all 31 bibliography entries (`ref-15`..`ref-104`, `ref-SI1`..`ref-SI10`) and all `from_source` / `described_from` links — exact match.

### Every datum in `data.yaml` checked against `pages/page-NNN.png`
- `r_drop_baseline` (p. 3): Verified `"1"` `um` against `page-003.png` (`R_drop = 1 μm`).
- `vasp_tet_conc_baseline` (p. 3): Verified `"0.40"` `uM` (trailing zero kept) against `page-003.png`.
- `n_seed_filaments` (p. 3): Verified `"30"` against `page-003.png`.
- `l_seed_initial` (p. 3): Verified `"0.1"` `um` against `page-003.png`.
- `l_seg_actin` (p. 17): Verified `"100"` `nm` against `page-017.png` (and `0.1 μm` on `page-003.png`, `100 nm` on `page-010.png`).
- `elongation_rate_baseline` (p. 3): Verified `"10.3"` `nm/s` against `page-003.png`, `page-017.png`, and `0.0103 μm/s` on `page-020.png`.
- `t_sim_baseline` (p. 3): Verified `"600"` `s` (`from_source: ref-15`) against `page-003.png` and `page-017.png`.
- `dt_langevin` (p. 17): Verified `"0.002"` `s` against `page-017.png` (and `2 ms` on `page-009.png`–`page-010.png`).
- `dt_movie_frame` (p. 37): Verified `"5"` `s` against `page-037.png`.
- `time_to_diameter_2R` (p. 3): Verified `"184.47"` `s` against `page-003.png`.
- `final_actin_conc_baseline` (p. 3): Verified `"27.67"` `uM` against `page-003.png`.
- `viscosity_droplet_printed` (p. 17): Verified literal `"0.5"` `pN/(s*um^2)` against `page-017.png` (`0.5 pN/s μm^2`, misprinted unit).
- `viscosity_ratio_water_printed` (p. 17): Verified `"50"` (`50x water`) against `page-017.png`.
- `boundary_stiffness_actin` (p. 17): Verified `"200"` `pN/um` against `page-017.png`.
- `boundary_stiffness_vasp` (p. 17): Verified `"100"` `pN/um` against `page-017.png`.
- `brownian_ratchet_force` (p. 18): Verified `"10"` `pN` (`ref-SI3`) against `page-018.png`.
- `actin_flexural_rigidity` (p. 18): Verified `"0.075"` `pN*um^2` (`ref-SI4`) against `page-018.png`.
- `actin_steric_radius` (p. 18): Verified `"3.5"` `nm` against `page-018.png`.
- `actin_steric_stiffness` (p. 18): Verified `"1.0"` `pN/um` against `page-018.png`.
- `kd_vasp_actin` (p. 3): Verified `"1.8"` `nM` (`ref-18`) against `page-003.png`.
- `vasp_disordered_pct` (p. 10): Verified `"about 60"` `%` (`ref-15`) against `page-010.png`.
- `vasp_polymerization_speedup` (p. 5): Fixed from `"3"` to `"≤ 3"` (`ref-18`) against `page-005.png` (see Section 1).
- `vasp_radius` (p. 18): Verified `"30"` `nm` against `page-018.png`, `page-019.png`, and `page-010.png`.
- `vasp_diffusion_rate` (p. 18): Verified `"10"` `um^2/s` against `page-018.png`.
- `vasp_binding_distance` (p. 18): Verified `"30"` `nm` against `page-018.png` and `page-010.png`.
- `vasp_valency` (p. 18): Verified `"4"` against `page-018.png`.
- `vasp_characteristic_unbinding_force` (p. 19): Verified `"10"` `pN` (`ref-SI5`) against `page-019.png`.
- `vasp_steric_stiffness` (p. 19): Verified `"1.0"` `pN/um` against `page-019.png`.
- `k_bind_sweep_range` (p. 3): Verified `"1 x 10^-4 to 10"` `1/s` against `page-003.png` and `page-018.png`.
- `k_unbind_sweep_range` (p. 3): Verified `"1 x 10^-4 to 10"` `1/s` against `page-003.png` and `page-018.png`.
- `k_bind_ring_fig4` (p. 5): Verified `"10"` `1/s` against `page-005.png`.
- `k_unbind_ring_fig4` (p. 5): Verified `"1"` `1/s` against `page-005.png`.
- `k_bind_fig5_text` (p. 5): Verified `"0.1"` `1/s` against `page-005.png`.
- `k_bind_fig5_caption` (p. 6): Verified `"1.0"` `1/s` against `page-006.png`.
- `k_unbind_fig5` (p. 6): Verified `"0.1"` `1/s` against `page-006.png`.
- `t_2pir_range` (p. 5): Verified `"10 to 600"` `s` against `page-005.png` and `page-006.png`.
- `vasp_tet_conc_range` (p. 18): Verified `"0.2 to 2.0"` `uM` against `page-018.png`.
- `k_bind_fig6` (p. 7): Verified `"1.0"` `1/s` against `page-007.png`.
- `k_unbind_fig6_range` (p. 7): Verified `"0.1 to 1.75"` `1/s` against `page-007.png` (including `1.5625` x-axis typo).
- `nucleator_molecules_table` (p. 20): Verified all 8 concentration/count pairs (`0..150 nM` -> `0..378`) against `page-020.png`.
- `nucleator_threshold_conc` (p. 5): Verified `"50"` `nM` against `page-005.png`.
- `k_bind_nucleator` (p. 20): Verified `"0.2"` `1/s` (`ref-SI6`) against `page-020.png`.
- `k_unbind_nucleator` (p. 20): Verified `"0.005"` `1/s` (`ref-SI7`) against `page-020.png`.
- `total_actin_conc_nucleator` (p. 20): Verified `"138.38"` `uM` (`ref-SI9`) against `page-020.png`.
- `radial_shell_thickness` (p. 27): Verified `"25"` `nm` against `page-027.png`.
- `deformation_z_step` (p. 10): Verified `"10"` `nm` against `page-010.png`.
- `deformation_relax_time` (p. 10): Verified `"1"` `s` against `page-010.png`.
- `aspect_ratio_max_approx` (p. 10): Verified `"~5.56"` against `page-010.png`.
- `fig7aii_aspect_ratio_vs_time` (p. 8): Verified all 6 `(time, a/b)` pairs (`604.502..1028.502 s` -> `1.03..5.56`) against `page-008.png` (and integer rounded labels on `page-035.png`).
- `fig7c_representative_aspect_ratios` (p. 8): Verified all 10 final aspect ratios (`1.75` [printed `1,75`], `3.95`, `4.34`, `4.96`, `3.40`, `4.14`, `1.97`, `1.97`, `2.16`, `2.06`) and exact decimal conversions of `{2/3..2}*pi` um index values against `page-008.png` (and nearest-integer purple highlights on `page-035.png`).
- `cap_total_volume_fraction` (p. 10): Verified `"12.5"` `%` against `page-010.png` and `page-034.png`.
- `cap_single_volume_fraction` (p. 34): Verified `"6.25"` `%` against `page-034.png`.
- `cap_conc_threshold_percent` (p. 34): Verified `"33.34"` `%` against `page-034.png`.
- `icosphere_subdivision_level` (p. 10): Verified `"4"` against `page-010.png` and `page-023.png`.
- `surface_monomer_cutoff` (p. 10): Verified `"100"` `nm` against `page-010.png`.
- `pca_var_3pc_main_fig3` (p. 4): Verified `"97.64"` `%` against `page-004.png`.
- `pca_var_suppfig4_caption` (p. 25): Verified `"95.83"` `%` against `page-025.png`.
- `pca_var_suppfig9_caption` (p. 30): Verified `"99.06"` `%` against `page-030.png`.
- `pca_var_suppfig11_caption` (p. 32): Verified `"96.86"` `%` against `page-032.png`.
- `pca_var_suppfig12_caption` (p. 33): Verified `"94.72"` `%` against `page-033.png`.
- `supptable3_ring_pct_probabilistic` (p. 19): Verified `"66.67"` `%` at `(10, 0.1)` against `page-019.png`.
- `supptable3_shell_sb_pct_probabilistic` (p. 19): Verified `"33.33"` `%` at `(10, 0.1)` against `page-019.png`.
- `lammps_sigma_vasp_crowder` (p. 16): Verified `"1"` against `page-016.png`.
- `lammps_sigma_actin` (p. 16): Verified `"0.1"` against `page-016.png` and `page-021.png`.
- `lammps_sigma_actin_crowder` (p. 16): Verified `"0.55"` against `page-016.png`.
- `lammps_eps_vasp_vasp` (p. 16): Verified `"1.66"` (`ref-SI1`) against `page-016.png`.
- `lammps_eps_actin_vasp` (p. 16): Verified `"1.2"` against `page-016.png`.
- `lammps_eps_repulsive` (p. 16): Verified `"1.0"` against `page-016.png`.
- `lammps_cutoff_attractive` (p. 16): Verified `"2.5"` against `page-016.png` and `page-020.png`.
- `lammps_cutoff_actin_actin` (p. 16): Verified `"1.5"` against `page-016.png`.
- `lammps_k_bond_actin` (p. 16): Verified `"30000"` against `page-016.png`.
- `lammps_r0_bond_actin` (p. 16): Verified `"0.1"` against `page-016.png`.
- `lammps_k_bend_phys` (p. 16): Verified `"10394"` `pN*nm` (`ref-SI2`) against `page-016.png`.
- `kbt_printed` (p. 16): Verified `"4.11"` `pN*nm` against `page-016.png`.
- `lammps_k_bend_thermal_calc` (p. 16): Verified `"~2529"` against `page-016.png`.
- `lammps_k_angle_chosen` (p. 16): Verified `"2500"` (`origin: assumed`) against `page-016.png`.
- `lammps_theta0_angle` (p. 16): Verified `"180"` `deg` against `page-016.png`.
- `lammps_binding_site_sep` (p. 16): Verified `"1.5"` against `page-016.png`.
- `lammps_binding_distance` (p. 16): Verified `"1.2"` against `page-016.png` and `page-021.png`.
- `lammps_binding_attempt_steps` (p. 17): Verified `"10000"` against `page-017.png` and `page-021.png`.
- `lammps_binding_probability` (p. 17): Verified `"0.5"` against `page-017.png` and `page-021.png`.
- `lammps_max_bonds_per_site` (p. 17): Verified `"1"` against `page-017.png`.
- `lammps_k_crosslink` (p. 17): Verified `"100"` against `page-017.png`.
- `lammps_r0_crosslink` (p. 17): Verified `"1.1"` against `page-017.png`.
- `lammps_box_length` (p. 20): Verified `"60"` against `page-020.png`.
- `lammps_crowder_count` (p. 20): Verified `"41250"` against `page-020.png`.
- `lammps_crowder_vol_frac` (p. 20): Verified `"10"` `%` against `page-020.png`.
- `lammps_vasp_count` (p. 20): Verified `"20625"` against `page-020.png`.
- `lammps_init_sphere_radius_printed` (p. 20): Verified `"60"` against `page-020.png`.
- `lammps_reduced_temp` (p. 21): Verified `"1.0"` against `page-021.png`.
- `lammps_dt_reduced` (p. 21): Verified `"0.0005"` against `page-021.png`.
- `lammps_equilibration_steps` (p. 21): Verified `"2 x 10^6"` against `page-021.png`.
- `lammps_total_steps` (p. 22): Verified `"14 x 10^7"` against `page-022.png`.
- `lammps_interface_density_threshold` (p. 22): Verified `"0.1"` against `page-022.png`.
- `exp_vasp_monomer_conc` (p. 11): Verified `"20"` `uM` against `page-011.png` and `page-010.png`.
- `exp_vasp_tetramer_conc` (p. 11): Verified `"5"` `uM` against `page-011.png`.
- `exp_vasp_labeled_pct` (p. 10): Verified `"10"` `%` against `page-010.png`.
- `exp_actin_conc_table` (p. 11): Verified all 3 ratios (`20, 10, 7.5` -> `1, 2, 2.67` `uM`) against `page-011.png`.
- `exp_tris_conc` (p. 10): Verified `"50"` `mM` against `page-010.png`.
- `exp_buffer_ph` (p. 10): Verified `"7.4"` against `page-010.png`.
- `exp_nacl_conc` (p. 10): Verified `"150"` `mM` against `page-010.png`.
- `exp_tcep_conc` (p. 10): Verified `"5"` `mM` against `page-010.png`.
- `exp_peg8k_wv` (p. 10): Verified `"3"` `%` against `page-010.png`.
- `exp_polymerization_time` (p. 11): Verified `"15"` `min` against `page-011.png`.
- `exp_staining_time` (p. 11): Verified `"15"` `min` against `page-011.png`.
- `exp_replicates` (p. 7): Verified `"3"` against `page-007.png` and `page-011.png`.
- `exp_wt_disk_pct_1_7p5` (p. 7): Verified `"70 ± 13"` `%` (`ref-15`) against `page-007.png`.
- `exp_gabmut_disk_pct_1_7p5` (p. 7): Verified `"< 10"` `%` against `page-007.png`.
- `exp_gabmut_ring_pct_1_7p5` (p. 7): Verified `"53 ± 6"` `%` against `page-007.png`.
- `exp_gabmut_shell_pct_1_7p5` (p. 7): Verified `"37 ± 4"` `%` against `page-007.png`.
- `exp_fig7ai_actin_conc` (p. 8): Verified `"2"` `uM` (`ref-15`) against `page-008.png`.
- `exp_fig7ai_vasp_conc` (p. 8): Verified `"10"` `uM` (`ref-15`) against `page-008.png`.

---

## 5. Verification of All Models vs Printed Equations & `symbolic()` Consistency

All 38 closed-form models were verified both against the page images (`pages/page-NNN.png`) and
programmatically for `run()` vs `symbolic()` agreement (`< 1e-6` relative error):

- `actin_concentration_from_ratio` (p. 11): Verified `[Actin] = [VASP-mon] / r_{V:A}` (`from_words`), units & `symbolic()` match.
- `actin_nucleation_ode` (p. 20): Verified `runnable: false` with `q_actin_nucleation_ode_unstated` (`ref-SI9`).
- `bell_slip_bond_unbinding` (pp. 3, 10, 19): Verified `k_unbind = k_unbind0 * exp(f / f_u)`, `tau_res = 1 / k_unbind`, units & `symbolic()` match.
- `bending_stiffness_thermal_units` (p. 16): Verified `k_bend_phys / kBT` against `page-016.png` equation, `symbolic()` match.
- `brownian_noise_from_mobility` (p. 10): Verified `sqrt(2 * mu * k_B T * dt)` with `1e-3` `pN*nm -> pN*um` conversion against `page-010.png`, `symbolic()` match.
- `brownian_noise_stddev` (p. 10): Verified `sqrt(2 * D * dt)` against `page-010.png` (`\sqrt{2 D^i dt}`), `symbolic()` match.
- `brownian_ratchet_growth_rate` (pp. 9, 18): Verified `k_grow * exp(-f / f_grow)` (`from_words`, `ref-SI3`), `symbolic()` match.
- `cytosim_actin_vasp_droplet` (p. 10): Verified Eq. (1) `dx^{ij}(t) = mu f^{ij}_tot(t) dt + dB_j(t)` against `page-010.png`; `runnable: false` with `q_cytosim_not_runnable`.
- `droplet_2d_energy_balance` (pp. 8–9): Verified `runnable: false` with `q_droplet_2d_energy_unstated` (`ref-15`).
- `droplet_characteristic_lengths` (pp. 3, 5, 24): Verified `2 R`, `pi R`, `2 pi R`, `symbolic()` match.
- `droplet_viscosity_ratio_water` (p. 17): Verified `eta_drop / eta_water`, `symbolic()` match.
- `einstein_diffusion_relation` (p. 10): Verified `D = mu * k_B T` with `1e-3` unit conversion against `page-010.png`, `symbolic()` match.
- `ellipse_2d_constant_area` (pp. 8–9): Verified `b = R^2 / a`, `(a/b)_2D = (L_max / (2 R))^2` (gives `pi^2 = 9.8696` matching Fig. 7E `page-008.png`; symbolic `"pi^2"` cannot be parsed as a decimal datum in `precision.py`), `symbolic()` match.
- `elongation_rate_from_growth_time` (p. 18): Verified `1000 * (L_target - L0) / T_grow` (`nm/s`), `symbolic()` match.
- `filament_elongation_length` (p. 3): Verified `L0 + k_grow * 1e-3 * t`, `symbolic()` match.
- `filament_growth_time` (p. 3): Verified `(L_target - L0) / (k_grow * 1e-3)`, `symbolic()` match.
- `filament_length_limited_aspect_ratio` (pp. 8, 10): Verified `(L_max / (2 * R_drop))^(3/2)` against Fig. 7B caption & p. 10 (`pi^(3/2) = 5.5683`), `symbolic()` match.
- `final_actin_concentration` (p. 3): Verified `(N_fil * L_fil * 1000 / dl) / (N_A * (4/3) * pi * R^3 * 1e-15) * 1e6`, `symbolic()` match.
- `harmonic_repulsion_force` (pp. 17–19): Verified `k_spring * max(0, delta) * 1e-3`, `symbolic()` match.
- `icosphere_subdivision_geometry` (pp. 10, 23): Verified recursive icosahedron midpoint projection onto sphere (`n=4` -> `5120` faces, median edge `75.5 nm`, matching Supp. Fig. 2B `page-023.png`).
- `lammps_crowder_volume_fraction` (p. 20): Verified literal printed sphere-radius form `N * (4/3)*pi*r^3 / L^3 * 100` (`79.99%`), `symbolic()` match.
- `lammps_crowder_volume_fraction_diameter` (p. 20): Verified `variant_of` diameter form `N * (pi/6)*d^3 / L^3 * 100` (`9.999%`), `symbolic()` match.
- `lammps_droplet_aspect_ratio` (p. 22): Verified `c / a` against Supp. Fig. 1F y-axis (`page-022.png`), `symbolic()` match.
- `lammps_droplet_asphericity` (p. 22): Verified `(c^2 / (a * b))^(1/3)` against Supp. Fig. 1E y-axis (`page-022.png`), `symbolic()` match.
- `lammps_harmonic_angle_potential` (p. 16): Verified LAMMPS convention `K_angle * (theta - theta0)^2` with `deg -> rad` port conversion (`q_lammps_harmonic_half_factor_convention`), `symbolic()` match.
- `lammps_harmonic_bond_potential` (pp. 16–17): Verified LAMMPS convention `K_bond * (r - r0)^2`, `symbolic()` match.
- `lammps_lennard_jones_potential` (p. 20): Verified printed `4 * eps * [(sigma/r)^12 - (sigma/r)^6]` for `r <= r_cut` and `r_wca = 2^(1/6)*sigma` against `page-020.png`, `symbolic()` match; shift/cutoff ambiguity logged in `q_lammps_lj_shift_and_actin_cutoff_ambiguity`.
- `lammps_lorentz_berthelot_sigma` (p. 16): Verified `(sigma_1 + sigma_2) / 2` against `page-016.png`, `symbolic()` match.
- `lammps_md_time_unit` (p. 21): Verified literal misprinted `sigma * (m / eps_VV)` (units `s^2/m`) against `page-021.png`, `symbolic()` match.
- `lammps_md_time_unit_sqrt` (p. 21): Verified `variant_of` corrected `sigma * sqrt(m / eps_VV)` (units `s`), `symbolic()` match.
- `lammps_radial_number_density` (p. 22): Added during audit (`N / (4 * pi * r^2 * dr)` from Supp. Fig. 1B y-axis `page-022.png`), `symbolic()` match.
- `lammps_reduced_temp_from_lj_eps` (pp. 16, 21): Added during audit (`1 / eps_VV_lj`), `symbolic()` match.
- `lammps_reduced_temperature` (p. 21): Verified `kBT / eps_VV` against `page-021.png`, `symbolic()` match.
- `lammps_vasp_actin_droplet` (pp. 20–21): Verified `runnable: false` with `q_lammps_not_runnable`.
- `langevin_overdamped_1d` (p. 10): Verified Euler-Maruyama 1D Eq. (1) solver (`dt = 0.002 s`, `D_eff = var / (2 * t_actual)` -> `9.991 um^2/s` vs `10 um^2/s`).
- `mobility_from_diffusion` (p. 10): Verified `D / (kBT * 1e-3)`, `symbolic()` match.
- `molecule_count_from_concentration` (p. 20): Verified `C * 1e-9 * N_A * (4/3)*pi*R^3 * 1e-15`, `symbolic()` match.
- `monomer_depletion_elongation_rate` (p. 20): Verified printed `k_grow * N_fil * (1 - [F-actin]/[T-actin])` against `page-020.png`, `symbolic()` match.
- `morphology_fraction_complement` (p. 7): Verified `100 - P_ring - P_shell`, `symbolic()` match.
- `persistence_length_from_rigidity` (pp. 2, 18): Verified `kappa / (kBT * 1e-3)` (`18.25 um`), `symbolic()` match.
- `prolate_spheroid_cap_geometry` (pp. 10, 34): Verified exact spherical/spheroidal double-cap volume fraction `0.5 * (1 - u)^2 * (2 + u)` (independent of aspect ratio, `u_cap = 0.6955` at `12.5%`), `f_single = f_caps / 2`, `C_thresh = C_bulk / 3` (`9.22 uM`, matching Supp. Fig. 13C `page-034.png`).
- `prolate_spheroid_volume_conserving` (pp. 10, 17): Verified `b = sqrt(R^3/a)`, `a/b = (a/R)^(3/2)`, prolate spheroid surface area `2*pi*b^2*(1 + (a/(b*e))*arcsin(e))`, `symbolic()` match.
- `quasistatic_droplet_deformation_schedule` (pp. 8, 10): Verified `a(t) = R_drop + 0.5 * floor((t - t0)/dt_relax) * dz`, `(a/R)^(3/2)` against Fig. 7Aii (`page-008.png`), Supp. Fig. 13B (`page-034.png`), and Supp. Fig. 14 (`page-035.png`). Note: `run_quasistatic_step428_single` in `runs.yaml` has a slightly misleading id (it evaluates at `t = 600 s`, step 0, to feed `run_prolate_spheroid_at_t600`), which is harmless and preserved.
- `stokes_einstein_diffusion` (pp. 17–18): Verified `reconstructed` Stokes-Einstein formula `kBT / (6 * pi * eta * r)` (`0.0145 um^2/s` at `0.5 Pa*s`, `0.145 um^2/s` at `0.05 Pa*s`), `symbolic()` match.
- `total_actin_addition_rate` (p. 5): Verified `k_grow * N_fil`, `symbolic()` match.
- `vasp_tetramer_concentration` (p. 11): Verified `C_vasp_mon / valency`, `symbolic()` match.

---

## 6. Verification of `added.yaml`, Comparisons, Refinements, and Relations

### `added.yaml` (all 14 entries)
- `actin_monomer_rise_added` (p. 3): Proper physical basis (`2.7 nm`) and real alternatives (`2.7027`, `2.75`).
- `viscosity_water_added` (p. 17): Proper physical basis (`0.001 Pa*s`) and real alternatives (`0.00089`, `0.001002`).
- `viscosity_droplet_unit_corrected_added` (p. 17): Proper basis (`0.5 Pa*s`) and decisive alternative (`0.05 Pa*s`, flips `cmp_viscosity_ratio_water` from `not reproduced` to `reproduced`).
- `force_probe_bell_added` (pp. 3, 18): Appropriate probe load (`5 pN`, alternatives `0`, `10`).
- `displacement_probe_spring_added` (pp. 17–19): Appropriate probe overlap (`10 nm`, alternatives `1`, `30`).
- `lammps_pair_distance_probe_added` (p. 20): Appropriate interior LJ probe (`1.5`, alternatives `1.12246`, `2.0`).
- `lammps_bond_stretch_probe_added` (p. 16): Appropriate bond stretch probe (`0.11`, alternatives `0.105`, `0.12`).
- `lammps_angle_probe_added` (p. 16): Appropriate angle bend probe (`170 deg`, alternatives `175`, `160`).
- `lammps_particle_mass_si_added` (p. 21): Proper ~160 kDa VASP tetramer estimate (`2.66e-22 kg`, alternatives monomer/dimer-of-tetramers `6.64e-23`, `5.31e-22`).
- `lammps_sigma_si_added` (p. 21): Proper diameter scale (`60 nm`, alternatives `30`, `10`).
- `lammps_eps_si_added` (pp. 16, 21): Audit disagreement logged in `q_audit_lammps_eps_si_circular_reading` (circular baseline `4.11 pN*nm`; alternatives `2.48`, `6.82` genuinely flip `cmp_lammps_reduced_temp`).
- `lammps_ellipsoid_a_probe_added` (p. 22): Appropriate minor-axis probe (`17.5`, alternatives `16.5`, `18.0`).
- `lammps_ellipsoid_c_probe_added` (p. 22): Appropriate major-axis probe (`20.5`, alternatives `19.0`, `22.0`).
- `lammps_shell_count_probe_added` (p. 22): Added during audit (`3000`, alternatives `2500`, `3500`).

### Comparisons (23 total)
- `cmp_elongation_rate_derived` (p. 3): `reproduced` (`10.3053` vs `10.3 nm/s`).
- `cmp_time_to_2R_printed_k` (p. 3): `reproduced` (`184.466` vs `184.47 s`).
- `cmp_time_to_2R_chained_k` (p. 3): `reproduced` (`184.371` vs `184.47 s`, via input-rounding band on `L0 = 0.1`).
- `cmp_time_to_2piR_printed_k` (p. 3): `reproduced` (`600.309` vs `600 s`).
- `cmp_final_actin_conc_exact_2pi` (p. 3): `reproduced` (`27.6757` vs `27.67 uM`).
- `cmp_final_actin_conc_chained_L600` (p. 3): `reproduced` (`27.6616` vs `27.67 uM`).
- `cmp_nucleator_molecule_counts` (p. 20): `reproduced` (all 8 points `0..378`).
- `cmp_exp_vasp_tetramer_conc` (p. 11): `reproduced` (`5.0` vs `5 uM`).
- `cmp_exp_actin_conc_table` (p. 11): `reproduced` (`1.0, 2.0, 2.667` vs `1, 2, 2.67 uM`).
- `cmp_lammps_k_bend_thermal` (p. 16): `no verdict` (`2528.95` vs `"~2529"`).
- `cmp_lammps_sigma_actin_crowder` (p. 16): `reproduced` (`0.55` vs `0.55`).
- `cmp_lammps_crowder_vol_frac_printed` (p. 20): `not reproduced` (`79.99%` vs `10%`, radius `1 sigma` error).
- `cmp_lammps_crowder_vol_frac_diameter` (p. 20): `reproduced` (`9.999%` vs `10%`, diameter `1 sigma` variant).
- `cmp_lammps_reduced_temp` (p. 21): `reproduced` (`1.0` vs `1.0`, circular baseline; flips to `not reproduced` under alternatives `2.48` and `6.82`).
- `cmp_lammps_reduced_temp_vs_table1` (pp. 16, 21): `not reproduced` (`0.6024` vs `1.0`, added during audit).
- `cmp_aspect_ratio_max_length_limited` (p. 10): `no verdict` (`5.5683` vs `"~5.56"`).
- `cmp_fig7aii_deformation_schedule` (p. 8): `reproduced` (all 6 snapshots `1.03..5.56`).
- `cmp_cap_single_volume_fraction` (p. 34): `reproduced` (`6.25` vs `6.25%`).
- `cmp_cap_conc_threshold_percent` (p. 34): `close` (`33.3333%` vs `33.34%`).
- `cmp_viscosity_ratio_water` (p. 17): `not reproduced` (`500` vs `50`, flips to `reproduced` at `0.05 Pa*s`).
- `cmp_langevin_D_consistency` (p. 10): `reproduced` (`9.991 ± 0.136` vs `10 um^2/s`).
- `cmp_stokes_einstein_vs_vasp_D` (pp. 17–18): `not reproduced` (`0.0145` vs `10 um^2/s`, `reconstructed` model).
- `cmp_gabmut_disk_complement` (p. 7): `no verdict` (`10.0%` vs `"< 10%"`).

### Refinements (3 total) & Relations (4 total)
- `ref_langevin_dt` (p. 10): Settled (`9.991 -> 10.211 -> 10.111 -> 10.066`, shift `0.75% < 5%`).
- `ref_langevin_n_particles` (p. 10): Settled (`9.991 -> 10.093 -> 9.954 -> 9.995`, shift `0.04% < 5%`).
- `ref_cap_bisection_tol` (p. 10): Settled (`0.695465 -> 0.695454 -> 0.695453 -> 0.695453`, shift `0.0017% < 0.1%`).
- `rel_spheroid_to_length_limited` (p. 8): Verified `same_as` by `boxgraph relations`.
- `rel_monomer_depletion_to_constant_addition` (p. 20): Verified `reduces_to` (`4.9e-9`) by `boxgraph relations`.
- `rel_brownian_noise_einstein` (p. 10): Verified `same_as` by `boxgraph relations`.
- `rel_lammps_reduced_temp_lj` (p. 21): Verified `same_as` by `boxgraph relations` (added during audit).

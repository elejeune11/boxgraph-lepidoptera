# Audit Log: `nastkolb2022`

Independent second-agent audit of `corpora/actin_bundling/records/nastkolb2022/` against `corpus/actin_bundling/2022_NastKolb.pdf` (PDF pages 1–18, including Supplemental Materials pp. 13–18), verified directly against rendered page images (`pages/page-001.png`–`pages/page-018.png`) and extracted text (`text/page-001.txt`–`text/page-018.txt`).

---

## Fixes Made

- `actin_conc_standard` (p. 11): fixed quote from `"Actin was added to a final concentration of 2 µM, of which 12.5% was"` (spanned lines 11–12 across column 2 in `text/page-011.txt`) to single-column-line verbatim quote `"was added to a final concentration of 2 µM, of which 12.5% was"`.
- `actin_labeled_fraction_atto488` (p. 11): fixed quote from `"Actin was added to a final concentration of 2 µM, of which 12.5% was"` (spanned lines 11–12 across column 2 in `text/page-011.txt`) to single-column-line verbatim quote `"was added to a final concentration of 2 µM, of which 12.5% was"`.
- `q_unstated_vasp_surface_density` (p. 11): fixed quote from `"Unbound protein was then removed by washing with 20 times the chamber"` (spanned lines 4–5 across column 2 in `text/page-011.txt`) to single-column-line verbatim quote `"protein was then removed by washing with 20 times the chamber"`.
- `actin_conc_high_fold_change`, `actin_concentration_fold_change`, `run_actin_conc_fold_change`, `cmp_actin_conc_fold_change_repro` (p. 2): added missing printed derived datum (`printed: "2"`, from `"centration of 4 µM, which is double the amount compared with the"`), `from_words: true` processing model (`R_actin = c_actin_high / c_actin_std`), run (`4 uM / 2 uM = 2.0`), and `reproduction` comparison (resolving the previously unused `concentration_fold_change` quantity in `quantities.yaml`; verdict: `reproduced`).

---

## Disagreements & New Open Questions Logged (First-Reader Judgements Preserved)

- `q_chimeric_vasp_construct_and_ml_substrate` (p. 10, `kind: reading`): first reading omits two experimental-condition details from `data.yaml`/`quantities.yaml` — (1) "VASP" is chimeric human VASP carrying the high-affinity GAB domain of *Dictyostelium discoideum* VASP (Breitsprecher et al., 2011), not wild-type human VASP; (2) SLBs were formed on piranha-etched glass coverslips while MLs were formed on uncoated ibidi polymer µ-slides from 30% isopropanol.
- `q_frap_amplitude_vs_mobile_fraction_reading` (p. 11 / p. 7, `kind: reading`): disagreement on first reader's FRAP parameterization and metadata — (1) `frap_recovery_fit` / `frap_recovery_by_half_time` bind fit amplitude `a` to `vasp_mobile_fraction` and evaluate `run_frap_fit_ml_end` with `a = 0.61` (`y(1.3 min) = 0.546`), whereas Fig. 4C (p. 7) visually shows `VASP ML prePol` reaching `~0.61–0.62` at `t = 1.3 min` (consistent with asymptote `a ~ 0.69` at `t_1/2 = 0.399 min`, where `0.69*(1 - 2^(-1.3/0.399)) = 0.618`), suggesting printed `0.61` may be the finite-window endpoint rather than `a`; (2) `frap_t_half_vasp_slb` and `frap_t_half_vasp_ml` use `origin: measured` and `polymerization_state: prePol`, whereas p. 11 derives `t_1/2 = tau*ln(2)` from a curve fit (`fitted`/`derived` alternative) and p. 6 does not explicitly specify `prePol`; (3) fitting a single exponential to a 2D diffusive FRAP curve can induce a positive baseline offset `c` exceeding the `±0.02` alternatives in `frap_baseline_offset_added`.
- `q_max_bundle_width_3um_vs_fig3f_whisker` (p. 6 / p. 5, `kind: unexplained_difference`): text on p. 6 claims VASPΔTetra exhibits bundle widths `"up to 3.2 µm compared with 3 µm for tetrameric VASP (Figure 3F)"`, whereas in Fig. 3F (p. 5, whiskers = 1–99th percentile) the upper whisker for tetrameric VASP sits at `~2.5 µm` (matching Fig. 2E) while VASPΔTetra sits at `~3.2 µm`; only the Fig. 2D histogram tail for 0.25 µM tetrameric VASP reaches `3.0 µm`.
- `q_intensity_extrapolation_equation_reading` (p. 11, `kind: reading`): `bundle_filament_intensity_extrapolation` lists `equation: "N_f = I_bundle / I_single"`, whereas p. 11 describes the method only in words (`"integrating the intensity over the diameter of filament bundles and used the result to extrapolate the number of filaments for large bundles"`) without specifying whether extrapolation was a single-filament ratio or a regression across thin bundles.

---

## Verification of `data.yaml` Against Page Images (61 Entries)

- `bundle_width_median_vasp_025` (p. 6): verified `0.942` `um`, `measured`, single-line quote, conditions match `page-006.png`.
- `bundle_width_median_vasp_dtetra_025` (p. 6): verified `0.584` `um`, `measured`, single-line quote, conditions match `page-006.png`.
- `bundle_width_max_vasp_025` (p. 6): verified `3` `um`, `measured`, single-line quote matches `page-006.png` (also matches p. 2 `"up to 3 µm"`; Fig. 3F whisker discrepancy logged in `q_max_bundle_width_3um_vs_fig3f_whisker`).
- `bundle_width_max_vasp_dtetra_025` (p. 6): verified `3.2` `um`, `measured`, single-line quote matches `page-006.png`.
- `bundle_diffraction_limited_window` (p. 2): verified `0 to 80` `s`, `measured`, single-line quote matches `page-002.png`.
- `filament_count_1um_intensity_measured` (p. 11): verified `45 ± 12` `dimensionless`, `measured`, single-line quote matches `page-011.png`.
- `vasp_diameter_range` (p. 11): verified `28 to 70` `nm`, `reported`, single-line quote matches `page-011.png`.
- `vasp_diameter_min` (p. 11): verified `28` `nm`, `reported`, single-line quote matches `page-011.png`.
- `vasp_diameter_max` (p. 11): verified `70` `nm`, `reported`, single-line quote matches `page-011.png`.
- `bundle_width_ref_1um` (p. 11): verified `1` `um`, `reported`, single-line quote matches `page-011.png`.
- `filament_count_1um_expected_min` (p. 11): verified `14` `dimensionless`, `predicted`, single-line quote matches `page-011.png`.
- `filament_count_1um_expected_max` (p. 11): verified `40` `dimensionless`, `predicted`, single-line quote matches `page-011.png`.
- `frap_mobile_fraction_slb_prepol` (p. 6): verified `0.97` `dimensionless`, `measured`, single-line quote matches `page-006.png`.
- `frap_mobile_fraction_slb_pol30` (p. 6): verified `0.65` `dimensionless`, `measured`, single-line quote matches `page-006.png`.
- `frap_mobile_fraction_ml_prepol` (p. 8): verified `0.61` `dimensionless`, `measured`, single-line quote matches `page-008.png`.
- `frap_mobile_fraction_ml_pol30` (p. 8): verified `0.57` `dimensionless`, `measured`, single-line quote matches `page-008.png`.
- `frap_t_half_vasp_slb` (p. 6): verified `0.102` `min`, single-line quote matches `page-006.png` (origin/`prePol` reading note logged in `q_frap_amplitude_vs_mobile_fraction_reading`).
- `frap_t_half_vasp_ml` (p. 6): verified `0.399` `min`, single-line quote matches `page-006.png` (origin/`prePol` reading note logged in `q_frap_amplitude_vs_mobile_fraction_reading`).
- `vasp_mobility_reduction_ml_vs_slb` (p. 6): verified `4` `dimensionless`, `derived`, single-line quote matches `page-006.png`.
- `frap_post_pol_time` (p. 7): verified `30` `min`, `reported`, single-line quote matches `page-007.png` Fig. 4C caption.
- `vasp_llps_threshold_main` (p. 6): verified `≥ 0.5` `uM`, `measured`, single-line quote matches `page-006.png`.
- `vasp_llps_range_si` (p. 14): verified `0.25 to 2` `uM`, `measured`, single-line quote matches `page-014.png` Supplemental Fig. 4G caption.
- `lipid_relaxation_fraction_60min` (p. 8): verified `~0.5` `dimensionless` (printed in words `"around half"`), `measured`, single-line quote matches `page-008.png`.
- `mdia1_tirf_decline_onset_time` (p. 6): verified `6` `min`, `measured`, single-line quote matches `page-006.png`.
- `network_overlay_time` (p. 2): verified `5` `min`, `reported`, single-line quote matches `page-002.png`.
- `actin_exchange_time` (p. 4): verified `3` `min`, `reported`, single-line quote matches `page-004.png` Fig. 2F caption.
- `actin_conc_standard` (p. 11): verified `2` `uM`, `reported` against `page-011.png`; cross-line quote fixed (see Fixes Made).
- `actin_conc_low_05` (p. 6): verified `0.5` `uM`, `reported`, single-line quote matches `page-006.png`.
- `actin_conc_mid_1` (p. 14): verified `1` `uM`, `reported`, single-line quote matches `page-014.png` Supplemental Fig. 3 caption.
- `actin_conc_high_4` (p. 2): verified `4` `uM`, `reported`, single-line quote matches `page-002.png`.
- `actin_conc_high_fold_change` (p. 2): added by audit (`2` `dimensionless`, `derived`, from `"double"` on `page-002.png`, single-line quote verified).
- `actin_conc_exchange_atto647` (p. 4): verified `0.5` `uM`, `reported`, single-line quote matches `page-004.png` Fig. 2F caption (Atto647N confirmed on p. 5 & p. 10).
- `actin_labeled_fraction_atto488` (p. 11): verified `12.5` `%`, `reported` against `page-011.png`; cross-line quote fixed (see Fixes Made).
- `cp_conc_standard` (p. 11): verified `50` `nM`, `reported`, single-line quote matches `page-011.png`.
- `profilin_conc_standard` (p. 11): verified `10` `uM`, `reported`, single-line quote matches `page-011.png`.
- `kcl_conc_kmei` (p. 11): verified `50` `mM`, `reported`, single-line quote matches `page-011.png`.
- `mgcl2_conc_kmei` (p. 11): verified `1` `mM`, `reported`, single-line quote matches `page-011.png`.
- `atp_conc_assay` (p. 11): verified `1` `mM`, `reported`, single-line quote matches `page-011.png`.
- `kmei_buffer_ph` (p. 11): verified `7.0` `dimensionless` (trailing zero preserved), `reported`, single-line quote matches `page-011.png`.
- `vasp_conc_005` (p. 2): verified `50` `nM`, `reported`, single-line quote matches `page-002.png`.
- `vasp_conc_0125` (p. 2): verified `0.125` `uM`, `reported`, single-line quote matches `page-002.png`.
- `vasp_conc_025` (p. 2): verified `0.25` `uM`, `reported`, single-line quote matches `page-002.png`.
- `vasp_conc_050` (p. 2): verified `0.5` `uM`, `reported`, single-line quote matches `page-002.png` (SI p. 14 misprint `0.5 M` already logged in `q_supplemental_fig2d_typo_05M`).
- `vasp_conc_100` (p. 2): verified `1` `uM`, `reported`, single-line quote matches `page-002.png`.
- `vasp_conc_bulk_050` (p. 14): verified `0.5` `uM`, `reported`, single-line quote matches `page-014.png` Supplemental Fig. 1D caption.
- `vasp_conc_beads_100` (p. 14): verified `1` `uM`, `reported`, single-line quote matches `page-014.png` Supplemental Fig. 1E caption.
- `vasp_dtetra_conc_025` (p. 5): verified `0.25` `uM`, `reported`, single-line quote matches `page-005.png` Fig. 3A caption.
- `alpha_actinin_conc_050` (p. 6): verified `0.5` `uM`, `reported`, single-line quote matches `page-006.png`.
- `alpha_actinin_conc_005` (p. 5): verified `0.05` `uM`, `reported`, single-line quote matches `page-005.png` Fig. 3B caption.
- `mdia1_conc_025` (p. 5): verified `0.25` `uM`, `reported`, single-line quote matches `page-005.png` Fig. 3C caption.
- `alpha_actinin_conc_combined_025` (p. 5): verified `0.25` `uM`, `reported`, single-line quote matches `page-005.png` Fig. 3D caption.
- `lipid_frac_egg_pc` (p. 10): verified `87.45` `%`, `reported`, single-line quote matches `page-010.png`.
- `lipid_frac_dogs_ni_nta` (p. 10): verified `10` `%`, `reported`, single-line quote matches `page-010.png`.
- `lipid_frac_dope_peg2000` (p. 10): verified `2.5` `%`, `reported`, single-line quote matches `page-010.png`.
- `lipid_frac_fluorescent_tracer` (p. 10): verified `0.05` `%`, `reported`, single-line quote matches `page-010.png` (four lipid fractions sum to exactly `100.00%`).
- `suv_extrusion_pore_size` (p. 10): verified `100` `nm`, `reported`, single-line quote matches `page-010.png`.
- `flow_chamber_vol` (p. 10): verified `50` `uL`, `reported`, single-line quote matches `page-010.png`.
- `tirf_objective_na` (p. 11): verified `1.47` `dimensionless`, `reported`, single-line quote matches `page-011.png`.
- `frap_bleach_circle_diameter` (p. 11): verified `4.96` `um`, `reported`, single-line quote matches `page-011.png`.
- `n_bundles_intensity_avg` (p. 11): verified `5` `dimensionless`, `reported`, single-line quote matches `page-011.png`.
- `n_frames_lipid_coloc_avg` (p. 11): verified `10` `dimensionless`, `reported`, single-line quote matches `page-011.png`.

---

## Verification of Models Against Printed Equations & Text (17 Models)

- `frap_recovery_fit` (p. 11): verified `y = a*(1 - exp(-t/tau)) + c` against `page-011.png` (`y = a*(1 – e^(-t/τ)) + c`); sign, symbols, `model.py` (`run` and `symbolic`) match exactly (amplitude-vs-mobile-fraction reading note in `q_frap_amplitude_vs_mobile_fraction_reading`).
- `frap_half_time` (p. 11): verified `t_1/2 = tau * ln(2)` against `page-011.png` (`t_{1/2} = τ*ln 2`); `model.py` (`run` and `symbolic`) exact.
- `frap_tau_from_half_time` (p. 11): verified algebraic inverse `tau = t_1/2 / ln(2)` (`from_words: true`, `built_from: [frap_half_time]`); `model.py` exact.
- `frap_recovery_by_half_time` (p. 11): verified substitution `y = a*(1 - exp(-t*ln(2)/t_1/2)) + c` (`from_words: true`, `built_from: [frap_recovery_fit, frap_half_time]`); `model.py` exact and verified `same_as` `frap_recovery_fit`.
- `frap_mobile_fraction_plateau` (p. 11): verified asymptotic limit `y_inf = a + c` (`from_words: true`, `built_from: [frap_recovery_fit]`); `model.py` exact and verified `reduces_to` as `t -> oo`.
- `frap_mobility_ratio` (p. 6): verified `R_mob = t_1/2,ML / t_1/2,SLB` (`from_words: true`); `model.py` exact; reproduces printed `"4"` (`0.399 / 0.102 = 3.912`).
- `bundle_filament_count_2d_spacing` (p. 11): verified `N_f = W / d_spacing` (`from_words: true`, `1000 * W / d_spacing` for `W` in `um` and `d_spacing` in `nm`); reproduces printed `14` (`14.29`) and honestly yields `not reproduced` against printed `40` (`35.71`).
- `bundle_effective_spacing` (p. 11): verified inverse `d_spacing = W / N_f` (`from_words: true`, `1000 * W / n_f`); `model.py` exact.
- `bundle_width_vasp_scaling` (p. 2): verified `W ~ c_VASP` (`from_words: true`, `prefactor: unknown`, unit prefactor `1` in `model.py`); matches `record_guide.md` prefactor-unknown rule.
- `cp_to_actin_ratio` (p. 2): verified `R_CP = c_CP / c_actin` (`from_words: true`, `1e-3 * c_cp / c_actin` for `c_cp` in `nM` and `c_actin` in `uM`); `model.py` exact.
- `actin_concentration_fold_change` (p. 2): added by audit (`R_actin = c_actin_high / c_actin_std`, `from_words: true`, with `symbolic()`); reproduces printed `"double"` (`4 / 2 = 2.0`).
- `bilayer_actin_decoration` (p. 11): verified `runnable: false` with `f_dec = A_bundle / A_total` matching p. 11 text and linked to `q_nonrunnable_image_processing_models`.
- `bundle_filament_intensity_extrapolation` (p. 11): verified `runnable: false` and linked to `q_nonrunnable_intensity_extrapolation`; reading note on inferred equation form logged in `q_intensity_extrapolation_equation_reading`.
- `bundle_intensity_ratio` (p. 11): verified `runnable: false` with `R_I = I_under / I_next` matching p. 11 text and linked to `q_nonrunnable_image_processing_models`.
- `bundle_width_lognormal_fit` (p. 4): verified `runnable: false` matching Fig. 2D / S2I / S4C captions and linked to `q_nonrunnable_lognormal_fit`.
- `frap_bleach_normalization` (p. 11): verified `runnable: false`, `described_from: ref-siggia2000`, matching p. 11 text and linked to `q_nonrunnable_image_processing_models`.
- `vasp_actin_spearman_colocalization` (p. 11): verified `runnable: false` matching p. 11 text and linked to `q_nonrunnable_image_processing_models`.

---

## Verification of `added.yaml`, Runs, Comparisons, Relations, Metadata & References

- `frap_baseline_offset_added` (p. 11 / p. 7): verified basis (`y(0) = 0` after Siggia et al. 2000 first-bleach-frame normalization, matching Fig. 4C) and alternatives `["-0.02", "0.02"]` (caveat on possible larger positive offset from exponential fit to 2D diffusion logged in `q_frap_amplitude_vs_mobile_fraction_reading`).
- `frap_observation_end_added` (p. 7): verified `1.3 min` and alternatives `["1.2", "1.4"]` directly against Fig. 4C x-axis on `page-007.png` (axis extends slightly past `1.2` to `~1.32 min`).
- `runs.yaml` & `comparisons.yaml` (pp. 2, 6, 8, 11): verified all 22 runs (21 original + 1 added) and 9 comparisons (8 original + 1 added); single validation warning (`run_frap_mobility_ratio: inputs disagree on substrate_type: ML | SLB`) confirmed benign and inherent to the cross-substrate ratio.
- `relations.yaml` (p. 11): verified all 3 stated relations (`rel_frap_fit_same_as_half_time`, `rel_frap_fit_reduces_to_plateau`, `rel_frap_half_time_fit_reduces_to_plateau`) pass `boxgraph relations` cleanly (`verified`).
- `source.yaml` (pp. 1–18): verified title, authors, year (`2022`), venue (`Molecular Biology of the Cell 33, ar91, 1–12`), DOI (`10.1091/mbc.E21-11-0577`), page map (`1–18`), and figure inventory (including broken y-axis ranges on Figs. 2E, 3F, S4F and tri-axis ranges on Fig. 5C) against `pages/page-001.png`–`pages/page-018.png`.
- `quantities.yaml`: verified dimensions and definitions of all 55 declared quantities (`bundle_pore_size` and `ni_nta_presence` reserved for upcoming figure digitization; `concentration_fold_change` now wired up).
- `references.yaml` (pp. 11–12): verified all 66/66 printed references present in exact order with `ref-claessens2008` linked to `source: claessens2008`.
- `questions.yaml` (pp. 2–14): verified all 15 original open questions (including Soumpasis diffusivity arithmetic `D_SLB ~ 0.225 um^2/s`, `D_ML ~ 0.058 um^2/s`; Rayleigh `0.216 um` / Abbe `0.177 um`; all four panel misreferences and SI `0.5 M` typo) against page images; fixed 1 cross-line quote (`q_unstated_vasp_surface_density`) and added 4 new audit questions.

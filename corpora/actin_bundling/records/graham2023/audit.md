# Independent Audit Log — `graham2023`

Second, independent reading of Graham et al., *Nature Physics* **19**, 574–585 (2023) (`corpus/actin_bundling/2023_Graham.pdf`, PDF pages 1–44: Main text + Methods + Extended Data Figs. 1–9 + Reporting Summary + Supplementary Information `MOESM1_ESM.pdf`) against the BoxGraph record `corpora/actin_bundling/records/graham2023/`. Every datum in `data.yaml` and every equation in `models/` was verified directly against the rendered page images (`pages/page-NNN.png`) and page text (`text/page-NNN.txt`). Neither `graham_code` nor `graham_source_data` was consulted, and figure digitization was left to the downstream stage.

---

## 1. Genuine Transcription Fixes Made

- `monomer_addition_rate_text` (p. 9 / printed p. 582): Fixed `printed` in `data.yaml` from `"8 to 9"` to `"about 8 to 9"`. The page image `pages/page-009.png` literally prints `"about 8 to 9 monomers per second"` (also inside the datum's own quote). Dropping the printed `"about"` qualifier violates BoxGraph's approximate-value rule (cf. `vasp_disordered_frac: "about 60"` and `vasp_droplet_fusion_time: "about 1"`) and manufactured hard-interval `reproduced` verdicts for `cmp_monomer_rate_from_mean_repro` and `cmp_monomer_rate_from_fit_repro` instead of `no verdict`.
- `cmp_ring_rod_transition_ar_repro` (p. 8 / p. 581; p. 44 / SI-17): Fixed an inaccurate hardcoded computed digit in the prose `note` of `comparisons.yaml` (`"AR* = 9.5364"`), replacing it with a non-numerical description; `ring_rod_critical_aspect_ratio` actually computes `9.53620975` (exact elliptic integral root of $P(a/b) = 4\pi R$, while the 7-term Gauss-Kummer series gives `9.536638`).

---

## 2. New Comparisons Added

- `cmp_actin_lp_table1_vs_text` (p. 3 / p. 576; p. 39 / SI-12): Added `consistency` comparison verifying `actin_lp_table1` (`17.73 um`) against `actin_lp_range_text` (`10 to 20 um`).
- `cmp_droplet_diam_median_vs_range` (p. 3 / p. 576; p. 19): Added `consistency` comparison verifying `droplet_diam_median_20uM` (`2.86 um`, Extended Data Fig. 4e) against `droplet_diam_range_obs` (`1 to 7 um`, main text).
- `cmp_droplet_diam_max_vs_range` (p. 3 / p. 576; p. 19): Added `consistency` comparison checking `droplet_diam_max_20uM` (`7.22 um`, Extended Data Fig. 4e) against `droplet_diam_range_obs` (`1 to 7 um`, main text; rounds `0.6–7.22 um` to `1–7 um`, verdict `close`).

---

## 3. Disagreements (`kind: reading`) and New Open Questions Added

- `q_audit_slsqp_direct_init_vs_extension_protocol` (`kind: reading`, p. 7 / p. 580; p. 41 / SI-14): Disagreement with first-reader implementation scope in `droplet_ring_minimization_printed` and `droplet_ring_minimization_corrected`. Both bypass SI Section 2.6's 21-cycle filament extension ($L = D \to 3.1D$) and iterative $15\text{ nm}$ $T_{\text{ring}}$ sweep, initializing directly on a circle of radius `r_init_frac * R` with `r_init_frac: 0.96` hardcoded in settings rather than exposed in `added.yaml`. Independent sensitivity checks at $D = 1\,\mu\text{m}, \gamma_I = 5.0\text{ pN}, T_{\text{ring}} = 300\text{ nm}$ show non-convexity: `r_init_frac = 0.90–0.96` converges to `AR = 1.2435` (`E = 24834 pN nm`, invariant across seeds `0, 1, 2, 42, 99`), whereas `r_init_frac = 0.85` stops at `AR = 1.211` and `r_init_frac = 0.98–0.99` gets trapped in a higher-energy local minimum (`AR = 1.81`, `E = 25052 pN nm`). Neither `r_init_frac` nor `maxiter` is refined in `refinements.yaml`.
- `q_audit_ring_bending_true_polar_curvature` (`kind: reading`, p. 4 / p. 577; pp. 42–43 / SI-15–16): Disagreement/refinement on `ring_bending_energy_eq1` and `ring_bending_energy_eq1_corrected`. While `ring_bending_energy_eq1_corrected` faithfully integrates SI Eq. S15 to get $\frac{1}{2}(1 + 6m^2 + m^4)$, expanding the exact polar curvature $\kappa(\theta) = |r^2 + 2r_\theta^2 - r r_{\theta\theta}|/(r^2 + r_\theta^2)^{3/2}$ with $ds = \sqrt{r^2 + r_\theta^2}\,d\theta$ at fixed contour length $L = 2\pi R$ yields $\frac{1}{2}(m^2 - 1)^2 = \frac{1}{2}(1 - 2m^2 + m^4)$, which properly vanishes for the $m = 1$ rigid translation mode (unlike both Eq. S21's $2+m^4$ and Eq. S15's $(1+6m^2+m^4)/2$). Also logs that Eq. S15 drops the minus sign in $r'(s)$ and that Eq. 1 / Eq. S21 sum over $m$ while `ring_bending_energy_eq1` evaluates a single mode $m$.
- `q_audit_scaling_exponent_vs_experiment` (`kind: unexplained_difference`, p. 5 / p. 578; p. 6 / p. 579; p. 7 / p. 580): Documents that all three spherical/AR < 1.1 experimental scaling-exponent validations (`cmp_exponent_pred_vs_fig3f_ar_lt11`: `1.5` vs `1.8037`; `cmp_exponent_pred_vs_fig4h_ar_lt11`: `1.5` vs `1.7195`; `cmp_exponent_pred_vs_fig3e_rings`: `1.5` vs `1.9682`) return `not reproduced` with no prior question in `questions.yaml`, and no SE/CI is printed on the empirical power-law exponents $b$.
- `q_audit_fig5d_ar_exceeds_pi_squared` (`kind: possible_error`, p. 8 / p. 581; p. 44 / SI-17): Fig. 5c caption and SI Section 2.11 state the maximum permitted aspect ratio is $(a/b)_{\max} = \pi^2 \approx 9.8696$ (where $2a = L = 2\pi R$), yet Fig. 5d plots configurations at `9.9` and `10.0` (and Fig. 5c extends to `10.0`), exceeding $\pi^2$.
- `q_audit_unstated_diffraction_limit_and_wetting_sweep_params` (`kind: unstated_input`, pp. 2, 9, 14–15, 21, 40): Logs unstated diffraction-limit subtraction value in Methods p. 15, unstated $T_{\text{ring}}$, $\gamma_I$, and initial geometry for the SI Section 2.4.2 / Extended Data Fig. 6 wetting sweep (straight $1.2\,\mu\text{m}$ chord in $1\,\mu\text{m}$ circle exposes $200\text{ nm}$ total vs `"100 nm"` / wetted fraction $\sim 0.917$), unstated SD vs SEM for FRAP/elongation `±` values, `"radius (or ring thickness)"` wording in SI Section 2.4, and `"Figs. 6 and 7"` typo on Methods p. 14.

---

## 4. Verification of `data.yaml` Against Page Images (50 Datums)

- `vasp_sequence_len` (p. 16, Extended Data Fig. 1; also Fig. 1a p. 2): Verified `"380"` dimensionless, quote verbatim. OK.
- `vasp_disordered_frac` (p. 2 / p. 575): Verified `"about 60"` percent, quote verbatim, refs 15–17 (`from_source: ref-15`). OK.
- `vasp_droplet_fusion_time` (p. 2 / p. 575): Verified `"about 1"` s, droplet buffer conditions match Fig. 1b caption, quote verbatim. OK.
- `frap_recovery_fig1f` (p. 2 / p. 575; Methods p. 14): Verified `"93.9 ± 4.89"` percent, 15 uM VASP, 3% Atto-488, 1.8 um ROI, quote verbatim. OK.
- `frap_recovery_text` (p. 2 / p. 575): Verified `"94 ± 5"` percent, quote verbatim. OK.
- `frap_t_half` (p. 2 / p. 575): Verified `"18 ± 3"` s, quote verbatim. OK.
- `droplet_diam_range_obs` (p. 3 / p. 576): Verified `"1 to 7"` um, 20 uM VASP droplet buffer (citing Supplementary Fig. S4 / Extended Data Fig. 4e), quote verified. OK.
- `droplet_diam_median_20uM` (p. 19, Extended Data Fig. 4e; p. 32, SI Fig. S4e): Verified `"2.86"` um, 20 uM VASP, quote verbatim. OK.
- `droplet_diam_max_20uM` (p. 19, Extended Data Fig. 4e; p. 32, SI Fig. S4e): Verified `"7.22"` um (wraps onto line 16 of right column on p. 19; single-line quote on line 15 avoids column interleaving, clarified in note). OK.
- `deformed_ar_threshold` (p. 3 / p. 576): Verified `"1.1"` dimensionless, quote verbatim. OK.
- `linear_droplet_ar_high` (p. 6 / p. 579): Verified `"> 15"` dimensionless, 10 uM VASP / 2 uM actin (1:5), quote verified. OK.
- `fit_fig3e_rings_a` (p. 5 / p. 578; Methods p. 14): Verified `"0.0476"`, 20 uM VASP phalloidin ring-thickness assay, quote verbatim. OK.
- `fit_fig3e_rings_b` (p. 5 / p. 578): Verified `"1.9682"`, quote verbatim. OK.
- `fit_fig3e_rings_r2` (p. 5 / p. 578): Verified `"0.7278"`, quote verbatim. OK.
- `fit_fig3e_discs_a` (p. 5 / p. 578): Verified `"0.0797"`, quote verbatim. OK.
- `fit_fig3e_discs_b` (p. 5 / p. 578): Verified `"1.8456"`, quote verbatim. OK.
- `fit_fig3e_discs_r2` (p. 5 / p. 578): Verified `"0.8429"`, quote verbatim. OK.
- `fit_fig3f_ar_lt11_a` (p. 5 / p. 578): Verified `"0.0629"`, quote verbatim. OK.
- `fit_fig3f_ar_lt11_b` (p. 5 / p. 578): Verified `"1.8037"`, quote verbatim. OK.
- `fit_fig3f_ar_lt11_r2` (p. 5 / p. 578): Verified `"0.7111"`, quote verbatim. OK.
- `exp_ring_thickness_exponent_text` (p. 6 / p. 579): Verified `"1.8"`, quote verbatim. OK.
- `fit_fig3f_ar_gt11_a` (p. 5 / p. 578): Verified `"0.0917"`, quote verbatim. OK.
- `fit_fig3f_ar_gt11_b` (p. 5 / p. 578): Verified `"1.7315"`, quote verbatim. OK.
- `fit_fig3f_ar_gt11_r2` (p. 5 / p. 578): Verified `"0.8979"`, quote verbatim. OK.
- `fit_fig4h_ar_lt11_a` (p. 7 / p. 580): Verified `"0.0721"`, quote verbatim. OK.
- `fit_fig4h_ar_lt11_b` (p. 7 / p. 580): Verified `"1.7195"`, quote verbatim. OK.
- `fit_fig4h_ar_lt11_r2` (p. 7 / p. 580): Verified `"0.6990"`, quote verbatim. OK.
- `lit_ring_thickness_exponent_range` (p. 6 / p. 579): Verified `"1.5 to 2"`, refs 27, 30, 31 (`from_source: ref-27`), quote verbatim. OK.
- `elongation_rate_fit_fig5h` (p. 8 / p. 581; Methods p. 14): Verified `"0.0227"` um/s, 10 uM VASP / 2 uM actin, no EDTA/EGTA, 1 mM ATP, quote verbatim. OK.
- `elongation_fit_r2_fig5h` (p. 8 / p. 581): Verified `"0.9965"`, quote verbatim. OK.
- `elongation_rate_mean_text` (p. 9 / p. 582): Verified `"23 ± 3"` nm/s, quote verbatim. OK.
- `monomer_addition_rate_text` (p. 9 / p. 582): **FIXED** from `"8 to 9"` to `"about 8 to 9"` 1/s (see Section 1).
- `actin_lp_range_text` (p. 3 / p. 576): Verified `"10 to 20"` um, quote verbatim. OK.
- `actin_lp_table1` (p. 39 / SI-12, Table 1): Verified `"17.73"` um, `from_source: ref-si-9`, quote verified. OK.
- `droplet_diam_sim_range_table1` (p. 39 / SI-12, Table 1): Verified `"1 to 4"` um, quote verified. OK.
- `droplet_diam_sim_range_main` (p. 6 / p. 579): Verified `"0.5 to 4"` um, quote verbatim. OK.
- `filament_spacing_vasp_delta` (p. 39 / SI-12, Table 1; p. 41 / SI-14): Verified `"15"` nm, `from_source: ref-si-8`, quote verbatim. OK.
- `filament_spacing_ena` (p. 39 / SI-12, Table 1; p. 41 / SI-14): Verified `"12"` nm (printed without `≈` in Table 1; `≈ 12 nm` cited to [13] on p. 41), quote verbatim. OK.
- `segment_length_lcyl` (p. 39 / SI-12, Table 1): Verified `"100"` nm, quote verified. OK.
- `single_filament_k_stretch` (p. 39 / SI-12, Table 1): Verified `"100"` pN/nm, `from_source: ref-si-10`, quote verified. OK.
- `gamma_avr_table1` (p. 39 / SI-12, Table 1): Verified `"1 x 10^-8"` pN/nm, quote verbatim. OK.
- `gamma_avd_table1` (p. 39 / SI-12, Table 1): Verified `"1 x 10^4"` pN/nm, quote verbatim. OK.
- `interface_steepness_ke` (p. 40 / SI-13, Section 2.4): Verified `"64"` dimensionless, quote verified. OK.
- `wetting_sweep_bundle_len` (p. 40 / SI-13, Section 2.4.2): Verified `"1.2"` um, quote verbatim. OK.
- `wetting_sweep_droplet_diam` (p. 40 / SI-13, Section 2.4.2): Verified `"1"` um, quote verbatim (preserving paper's typo `"cirular"`). OK.
- `wetting_sweep_exposed_len` (p. 40 / SI-13, Section 2.4.2): Verified `"100"` nm, quote verbatim. OK.
- `ring_extension_step_frac` (p. 41 / SI-14, Section 2.6): Verified `"0.1"` dimensionless, quote verbatim. OK.
- `ring_extension_cycles_val` (p. 41 / SI-14, Section 2.6): Verified `"21"` dimensionless, quote verbatim. OK.
- `slsqp_tol_val` (p. 41 / SI-14, Section 2.6): Verified `"10^-6"` dimensionless, quote verified. OK.
- `noise_amp_val` (p. 41 / SI-14, Section 2.6): Verified `"0.1"` nm, quote verbatim. OK.
- `ring_thickness_step_val` (p. 7 / p. 580, Fig. 4e caption): Verified `"15"` nm, quote verbatim. OK.
- `fig4f_droplet_diam` (p. 7 / p. 580, Fig. 4f title): Verified `"1"` um (discrepancy with caption `"initial radius 1 um"` documented in `q_fig4f_radius_vs_diameter`), quote verbatim. OK.
- `fig4f_gamma_i` (p. 7 / p. 580, Fig. 4f title): Verified `"5.0"` pN, quote verbatim. OK.
- `fig4g_gamma_i_low` (p. 7 / p. 580, Fig. 4g caption): Verified `"1.0"` pN, quote verbatim. OK.
- `fig4h_gamma_i_max` (p. 7 / p. 580, Fig. 4h legend; p. 22, Extended Data Fig. 7g): Verified `"15.0"` pN, quote verified. OK.
- `fig4f_tring_min` (p. 7 / p. 580, Fig. 4f caption): Verified `"5"` nm, quote verbatim. OK.
- `fig4f_tring_max` (p. 7 / p. 580, Fig. 4f caption): Verified `"0.48"` um, quote verbatim. OK.
- `fig5c_droplet_radius` (p. 44 / SI-17, Section 2.11): Verified `"1"` um, quote verbatim. OK.
- `fig5c_kb_printed` (p. 8 / p. 581, Fig. 5c caption): Verified `"1"` pN*nm^2 (mismatch with SI Section 2.11 documented in `q_fig5c_unit_mismatch`), quote verbatim. OK.
- `fig5c_gamma_i_printed` (p. 8 / p. 581, Fig. 5c caption): Verified `"1"` pN, quote verbatim. OK.
- `pred_ring_thickness_exponent` (p. 6 / p. 579): Verified `"1.5"` dimensionless, quote verbatim. OK.
- `derived_final_ring_length_ratio` (p. 41 / SI-14, Section 2.6): Verified `"3.1"` dimensionless, quote verbatim. OK.
- `pred_ring_rod_transition_ar` (p. 44 / SI-17, Section 2.11; p. 8 / p. 581, Fig. 5d caption): Verified `"> 9"` dimensionless, quote verified. OK.

---

## 5. Verification of `models/` Against Printed Equations on Page Images (31 Models)

- `bundle_bending_rigidity_bathe_eqS1` (p. 39 / SI-12, Eq. S1): Matches printed $k_b = k_{b,\text{fil}} N [1 + \chi^2(N-1)/(1 + c(q_j)(N+\sqrt{N})/\alpha)]$, `run()` and `symbolic()` agree. OK.
- `bundle_bending_rigidity_linear` (p. 41 / SI-14, Section 2.6; p. 39 Note 1): Matches printed $K_{b,\text{ring}} = K_{b,\text{fil}} \times N_{\text{fil}}$. OK.
- `circular_ring_contour_length` (p. 8 / p. 581; p. 43 / SI-16): Matches printed $L = 2\pi R$. OK.
- `discrete_bending_constant` (p. 39 / SI-12, Table 1): Matches printed $k_{\text{bend}} = N_{\text{fil}} L_p k_B T / L_{\text{cyl}}$ with unit conversion $1\,\mu\text{m} = 1000\text{ nm}$. OK.
- `discrete_stretching_constant` (p. 39 / SI-12, Table 1): Matches printed $k_{\text{stretch}} = N_{\text{fil}} k_{\text{stretch,fil}}$. OK.
- `droplet_ring_minimization_printed` (pp. 39–41 / SI-12–14, Eqs. S2–S7, S11–S12): Faithfully implements literal printed Eq. S2 (no rest length $L_{\text{cyl}}$, $+k_{\text{bend}} u_i \cdot u_{i+1}$) and literal printed Eqs. S4–S5 ($\Delta = \frac{1}{2}(\gamma_{\text{AVR}} - \gamma_{\text{AVD}})$, $\gamma = \Delta \tanh + \Delta + \gamma_{\text{AVD}}$) with verified analytical gradients; initialization/protocol scope disagreement recorded in `q_audit_slsqp_direct_init_vs_extension_protocol`. OK.
- `droplet_ring_minimization_corrected` (pp. 7, 20, 39–41): `variant_of: droplet_ring_minimization_printed`. Verified analytical gradients of stretching, bending ($k_{\text{bend}}(1 - u_i \cdot u_{i+1})$, exact continuum limit $\frac{1}{2}k_b\int\kappa^2 ds$), wetting, and 7-term perimeter; reproduces Fig. 4f–h circle-to-ellipse deformation (`AR` crosses `1.1` around `T_ring = 90–120 nm`, `AR = 1.244` at `300 nm`, `1.926` at `480 nm` for `D = 1 um, gamma_I = 5 pN`); initialization/protocol scope disagreement recorded in `q_audit_slsqp_direct_init_vs_extension_protocol`. OK.
- `droplet_surface_energy_eqS12` (p. 41 / SI-14, Eq. S12; pp. 6–7): Matches printed $\Delta E_{\text{surface}} = \gamma_I L$ (`pN * um -> 1000 pN*nm`). OK.
- `effective_spherical_diameter` (p. 15, Methods): Matches printed $D = (a b^2)^{1/3}$. OK.
- `ellipse_axes_constant_area` (p. 39 / SI-12, Table 1; p. 44 / SI-17): Matches printed $4ab = D_{\text{droplet}}^2$ ($ab = R^2$) with $a/b = \text{AR}$. OK.
- `ellipse_perimeter_integral_eqS10` (p. 41 / SI-14, Eqs. S8, S10): Matches printed $L = 4\int_0^{\pi/2}\sqrt{(a\sin\theta)^2 + (b\cos\theta)^2}\,d\theta$. OK.
- `ellipse_perimeter_series_eqS11` (p. 41 / SI-14, Eq. S11): Matches printed 7-term Gauss-Kummer expansion $\pi(a+b)\sum_{n=0}^{6}\binom{0.5}{n}^2 h^n$, $h=((a-b)/(a+b))^2$. OK.
- `elongation_monomer_rate` (p. 9 / p. 582): Matches `from_words: true` conversion $r_{\text{mon}} = v / \delta l$. OK.
- `filament_flexural_rigidity` (p. 39 / SI-12, Note 1; p. 43 / SI-16): Matches printed $k_{b,\text{fil}} = L_p k_B T$ (`um * pN*nm -> 1000 pN*nm^2`). OK.
- `linear_bundle_elongation` (p. 8 / p. 581, Fig. 5h caption): Matches printed linear fit $y = v t$. OK.
- `power_law_thickness_fit` (p. 5 / p. 578, Fig. 3e–f; p. 7 / p. 580, Fig. 4h): Matches printed $y = a x^b$ in micrometre units. OK.
- `ring_amplitude_scaling_eqS24` (p. 4 / p. 577; p. 43 / SI-16, Eq. S24): Matches printed $u_m \approx R^2/\sqrt{L L_p}$ with `prefactor: unknown`. OK.
- `ring_bending_energy_eq1` (p. 4 / p. 577, Eq. 1; p. 42 / SI-15, Eq. S21): Matches printed $\frac{k_b}{2}(L/R^2 + (L/R^4)(2+m^4)u_m^2)$ for single mode $m$ (`1/um = 1e-3/nm`); single-mode vs sum-over-$m$ noted in `q_audit_ring_bending_true_polar_curvature`. OK.
- `ring_bending_energy_eq1_corrected` (p. 42 / SI-15, Eqs. S15–S21): `variant_of: ring_bending_energy_eq1`. Independently verified via SymPy that exact integration of SI Eq. S15 yields $\frac{1}{2}(1 + 6m^2 + m^4)$; exact polar curvature refinement recorded in `q_audit_ring_bending_true_polar_curvature`. OK.
- `ring_confinement_bending_energy` (p. 43 / SI-16, Section 2.8): Matches printed $\frac{k_b}{2} L/R^2$. OK.
- `ring_equipartition_amplitude_eqS23` (p. 43 / SI-16, Eq. S23): Matches algebraic solution $u_m = R^2/\sqrt{(2+m^4) L L_p}$ of Eq. S23 with $f(m) = 2+m^4$ and $L_p = k_b/k_B T$. OK.
- `ring_excess_bending_energy_eqS22` (p. 43 / SI-16, Eq. S22): Matches printed $\Delta E_{\text{bend}} = \frac{k_b}{2}(L/R^4) u_m^2 (2+m^4)$. OK.
- `ring_extension_length_ratio` (p. 41 / SI-14, Section 2.6): Matches printed $L_{\text{final}}/D = 1 + n_{\text{cycles}} \cdot \text{step\_frac}$. OK.
- `ring_filament_count` (p. 39 / SI-12, Table 1; p. 7 / p. 580, Fig. 4a): Matches printed $N_{\text{fil}} = T_{\text{ring}} / \delta$. OK.
- `ring_rod_critical_aspect_ratio` (p. 44 / SI-17, Section 2.11): Solves printed $P(R\sqrt{\text{AR}}, R/\sqrt{\text{AR}}) = 2L$ via `brentq` (`AR* = 9.53621`). OK.
- `ring_rod_max_aspect_ratio` (p. 8 / p. 581, Fig. 5c; p. 44 / SI-17): Matches printed $(a/b)_{\max} = L^2/(4R^2) = \pi^2$ when $2a = L = 2\pi R, ab = R^2$. OK.
- `ring_rod_transition_energy` (p. 44 / SI-17, Section 2.11; p. 8 / p. 581, Fig. 5c–d): Verified exact ellipse curvature integral $\kappa^2 ds/dt = (ab)^2 / ((a\sin t)^2 + (b\cos t)^2)^{5/2}$, minimum-placement search over $s_0 \in [0, P)$, and quantitative match to Fig. 5c plotted curves when $k_b = 1\text{ pN }\mu\text{m}^2, \gamma_I = 1\text{ pN}, R = 1\,\mu\text{m}$ (`E_surf(1) = 6.283`, `E_bend(1) = 3.142`, peak `E_bend(9.536) = 19.77`, plunging to `0.66` at `AR = 10`). OK.
- `ring_thickness_scaling_eqS25` (p. 43 / SI-16, Eq. S25): Matches printed $u_m \approx R^{1.5}/\sqrt{L_p}$ with `prefactor: unknown`. OK.
- `ring_thickness_scaling_exponent` (p. 6 / p. 579; p. 43 / SI-16): Computes $d(\ln u_m)/d(\ln D) = 1.5$ from Eq. S24 with $L = 2\pi R$ and constant $L_p$. OK.
- `wetting_tension_profile_eqS5` (p. 40 / SI-13, Eqs. S4–S5; p. 20, Extended Data Fig. 5c): Matches printed $\Delta = \frac{1}{2}(\gamma_{\text{AVR}} - \gamma_{\text{AVD}})$, $\gamma(e) = \Delta\tanh(k_e(e-1)) + \Delta + \gamma_{\text{AVD}}$. OK.
- `wetting_tension_profile_eqS5_corrected` (p. 20, Extended Data Fig. 5c): `variant_of: wetting_tension_profile_eqS5`. Matches plotted curve rising from $\gamma_{\text{AVR}}$ inside ($e < 1$) to $\gamma_{\text{AVD}}$ outside ($e > 1$). OK.

---

## 6. Audit of `added.yaml` (17 Entries)

- `kbt_room_temp_added` (pp. 14, 39, 43): Basis (298 K) and alternatives (`4.045, 4.142, 4.280` pN nm spanning 293–310 K) are physical and well-justified. OK.
- `actin_monomer_step_added` (p. 9 / p. 582): Basis (`2.7 nm` F-actin axial rise) and alternatives (`2.75, 2.73, 2.8 nm`) are physical and real. OK.
- `test_droplet_radius_added` (pp. 4, 19, 42–43): Half of median diameter `2.86 um` with alternatives `0.5, 1.0, 2.0 um`. OK.
- `test_fourier_mode_added` (pp. 4, 42–43): Mode `m = 2` with alternatives `1, 3, 4`. OK.
- `test_fluctuation_amp_added` (pp. 4–5, 42–43): `0.35 um` with alternatives `0.1, 0.5, 0.7 um`. OK.
- `test_ring_thickness_300nm_added` (p. 7, Fig. 4f): `300 nm` with alternatives `150, 240, 450 nm`. OK.
- `test_aspect_ratio_ellipse_added` (pp. 7–8, 41): `1.5` with alternatives `1.1, 2.0, 5.0`. OK.
- `test_aspect_ratio_high_added` (pp. 8, 44): `9.5` with alternatives `9.0, 9.6, 9.8696`. OK.
- `test_ellipse_coord_inside_added` (pp. 20, 40): `0.5` with alternatives `0.0, 0.25, 0.75`. OK.
- `test_ellipse_coord_outside_added` (pp. 20, 40): `1.5` with alternatives `1.25, 2.0, 2.5`. OK.
- `test_major_axis_methods_added` & `test_minor_axis_methods_added` (p. 15): `3.5 um` & `2.5 um` with real alternatives. OK.
- `test_bathe_alpha_weak_added`, `test_bathe_chi_added`, `test_bathe_cqj_added` (p. 39 / SI-12, Eq. S1): Weak-coupling test values (`0.01`, `1.5`, `1.0`) with plausible alternatives spanning orders of magnitude. OK.
- `test_elongation_time_added` (p. 8 / p. 581, Fig. 5h): `100 s` midpoint with alternatives `50, 150, 200 s`. OK.
- `fig5c_kb_um_consistent_added` (pp. 8, 44): `1 pN*um^2` with alternatives `0.5, 2.0, 10.0` (clearly distinguished from literal caption `fig5c_kb_printed = 1 pN*nm^2`). OK.

---

## 7. Verification of `refinements.yaml`, `relations.yaml`, and Existing `questions.yaml`

- `ref_ellipse_perim_integral_nquad`, `ref_ellipse_perim_series_nterms`, `ref_ring_rod_ar_crit_nquad`, `ref_ring_rod_energy_narc_ebend`, `ref_slsqp_printed_ftol`, `ref_slsqp_corrected_ftol`, `ref_slsqp_corrected_nwet`, `ref_slsqp_corrected_nperim` (pp. 40–44): All 8 refinements converge (`settles: True`); omission of `r_init_frac` and `maxiter` from refinement/added inputs is noted in `q_audit_slsqp_direct_init_vs_extension_protocol`. OK.
- `rel_bathe_reduces_to_linear`, `rel_eq1_reduces_to_confinement`, `rel_eq1_corrected_reduces_to_confinement`, `rel_ellipse_integral_approx_series`, `rel_eqS23_scales_as_eqS24`, `rel_eqS25_scales_as_eqS24`, `rel_ring_rod_reduces_to_circular_confinement` (pp. 39–44): All 7 stated intra-record relations verified via `boxgraph relations --only graham2023` (`verified`). OK.
- Prior questions `q_eqS2_stretch_bend_misprint`, `q_eqS4_S5_wetting_sign_inverted`, `q_fig5c_unit_mismatch`, `q_eqS15_S21_derivation_algebra`, `q_lp_bundle_vs_filament_scaling`, `q_fig4f_radius_vs_diameter`, `q_fig3f_vs_fig4h_ar_lt11_fits`, `q_unstated_kbt_and_temperature`, `q_unstated_wetting_quadrature_N`, `q_extension_step_100nm_vs_01D`, `q_methods_diameter_vs_2d_area_and_discs`, `q_unstated_monomer_step_and_surface_tension_dimensions` (pp. 5–9, 15, 20, 39–44): All 12 verified against the page images and preserved intact. OK.

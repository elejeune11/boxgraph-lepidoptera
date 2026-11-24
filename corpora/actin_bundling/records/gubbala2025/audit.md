# Audit log for `gubbala2025`

Independent page-image audit of `corpora/actin_bundling/records/gubbala2025/` against all 20 pages of `corpus/actin_bundling/2025_Gubbala.pdf` (PDF pp. 1–9 main text, pp. 10–20 Electronic Supplementary Information). All numbers in `data.yaml`, all 25 models in `models/` (algebraically verified with SymPy for Eqs. 3–7, 9–10, S1–S34), all 20 entries in `added.yaml`, all `references.yaml` entries, and all `relations.yaml` entries were independently checked against the page images.

## Transcription errors corrected

- `rupture_duration_si` (PDF p. 12 / SI-3, SI Sec. 1.5): changed `printed` from `"20"` to `"≥ 20"` because the page image prints `"heated to 37°C for at least 20 minutes"`, which is a lower bound rather than an exact 20 min duration.

## Missing printed results, runs, and comparisons added

- `S_Lo_printed` (PDF p. 3 / 6245, Sec. 2.1, Eq. 4): added the paper's printed prediction that the nematic order vanishes (`"0"`) in the no-actin Lo phase (`phi = -1`), quoting `"field becomes zero, indicating an isotropic phase."`.
- `run_nematic_sat_exact_eq4_Lo` (PDF p. 3 / 6245, Eq. 4): added run of the variant model `nematic_saturation_from_eq4_exact` at `phi = -1`, `beta = 1` to pair with existing `run_nematic_sat_beta1_Lo`.
- `cmp_S_Lo_zero` and `cmp_S_Lo_zero_exact` (PDF p. 3 / 6245, Eq. 4): added reproduction comparisons of `S_sat_beta1_Lo` and `S_sat_exact_eq4_Lo` against `S_Lo_printed` (both reproduced: `0.0` vs. `"0"`).
- `cmp_beta_min_lsa_corr_vs_default` (PDF p. 16 / SI-7, SI Sec. 3.3 Eq. S28): added consistency comparison of `betamin_lsac_low_x` (from variant `linear_stability_s24_s28_corrected`) against `beta_default` (`1.0` vs. `"1"`), matching existing `cmp_beta_min_lsa_vs_default`.
- `cmp_nondim_E3` (PDF p. 15 / SI-6, SI Sec. 3.2): added consistency comparison of `E3_nondim` against `E3_low_actin` (`0.1` vs. `"0.1"`), completing the dimensionless-group checks alongside `cmp_nondim_gamma`, `cmp_nondim_E1`, and `cmp_nondim_chi`.
- `cmp_pde_low_Rratio_printed` and `cmp_pde_low_Rratio_variational` (PDF p. 6 / 6248, Fig. 3): added reproduction comparisons of `Rratio_pde_low_printed` and `Rratio_pde_low_var` against `R_ratio_low_sim_fig3_end` (`~1.73`, `no verdict` until digitized), symmetrically matching existing `cmp_pde_high_Rratio_fig3`.

## Additional open questions added

- `q_eqS8_cayley_hamilton_3d_in_2d` (PDF pp. 13–14 / SI-4–SI-5, SI Sec. 3.1 Eq. S8): `possible_error` — Eq. (S8) states the 3D Cayley–Hamilton identity and p. 14 applies it to a 2D traceless matrix with `det Q = -(1/2) Tr(Q^2)`, producing the false 2D equation `Q^3 = (1/2) Tr(Q^2)(Q - I)` (whose trace is `-Tr(Q^2) != 0` instead of `0`).
- `q_high_actin_pde_fast_regime_mismatch` (PDF pp. 5, 6, 18 / 6247, 6248, SI-9, Sec. 3 / Fig. 3 / Fig. S2): `unexplained_difference` — `run_pde_high_variational_stab` over `t = 200..1200` produces `alpha = 0.405` and `R(1200)/R(200) = 2.04` (faster than both `1/4` and `1/3`, rather than `< 0.25` and `~1.32`) because that window falls inside the Fig. S2 `"fast"` transient (diagnostic evaluation over `t = 1000..3000` yields `alpha = 0.156 < 0.25` and `R_ratio = 1.19`), and `run_pde_high_variational_stab` has 1 replicate without a `dt` refinement.
- `q_abstract_vs_results_alpha_bound` (PDF pp. 1, 2, 5, 6 / 6243–6248, Abstract & Sec. 3): `overstated_claim` — Abstract (`alpha < 1/4`) and Introduction (`alpha < 0.25 for varying densities`) claim strict `< 1/4` across actin stiffnesses/densities, whereas Sec. 3 and Fig. 3 report `alpha approx 0.25` (consistent with `1/4`) at low actin density and `alpha < 0.25` only at high actin density.
- `q_si_methods_active_copy_paste` (PDF p. 11 / SI-2, SI Secs. 1.2–1.3): `possible_error` — SI Sec. 1.2 is titled `"Actin and myosin preparation"` yet describes only actin, and SI Sec. 1.3 states the second GUV recipe targeted `"the lowest area fraction sampled amongst the active samples"`, contradicting the passive (no-activity) study described in the main text.
- `q_audit_si_sec4_imex_prose_and_dt_exponent` (PDF p. 16 / SI-7, SI Sec. 4 Eqs. S31–S32): `possible_error` — prose states `"All linear terms are treated implicitly"` even though Eqs. (S31)–(S32) treat linear `-phi` and `+Q_ij` explicitly, and the adaptive step update writes exponent `1/5` (4th-order) for a 1st-order forward Euler IMEX scheme.

## Reading disagreements / judgements preserved (logged as `reading` questions)

- `q_audit_nondim_s32_lambda_prose` (PDF p. 15 / SI-6, `nondim_parameters_s32`): SI Sec. 3.2 prose literally prints the dimensionally mixed formula `lambda = lambda_tilde + chi/2`, whereas `nondim_parameters_s32` implements `lambda = 1 + chi/2` from dividing Eq. (S17) by `l_c^2 = lambda_tilde`. Preserved as the first reader's documented choice rather than modifying the model.
- `q_audit_bounds_S1_symbolic_eval` (PDF p. 3 / 6245, `E3_over_E1_min_S1`, `E3_over_E1_max_S1`): Page 3 prints symbolic bounds `-1.5*S <= E_3/E_1 <= 3/S` rather than explicitly printing the evaluated numbers `-1.5` and `3` at `S = 1`. Preserved as the first reader's choice under `conditions: {nematic_order_magnitude: "1"}`.

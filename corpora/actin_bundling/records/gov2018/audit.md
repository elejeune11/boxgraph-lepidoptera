# Audit of `gov2018`

Independent second reading of N. S. Gov, *Phil. Trans. R. Soc. B* 373: 20170115 (2018) (`corpus/actin_bundling/2018_Gov.pdf`, PDF pages 1–13) against `corpora/actin_bundling/records/gov2018/`. Every page image (`pages/page-001.png`–`pages/page-013.png`) and extracted text file (`text/page-001.txt`–`text/page-013.txt`) was inspected directly alongside every YAML file and `model.py` in the record.

## 1. Data Audit (`data.yaml` vs. page images)

- `bending_modulus_kBT_range` (p. 3): verified on `pages/page-003.png` col. 2 (`κ/k_B T ~ 1–100`); digits, unit (`dimensionless`), page, and single-column quote exact.
- `curvature_radius_ratio_range` (p. 3): verified on `pages/page-003.png` col. 2 (`a/R ~ 0.1–10^-3`); recorded as `"~0.1 to 0.001"` because `boxgraph/precision.py` cannot parse scientific exponents inside approximate ranges (documented in datum note and `q_audit_molecular_scales_abs_convention`).
- `continuum_validity_length_est` (p. 3): verified on `pages/page-003.png` col. 2 (`approximately 10 nm`); value `"≈ 10"`, unit `nm`, `from_source: ref-42` (noting `[42–45]`), and quote exact.
- `passive_active_force_zero` (p. 3): verified on `pages/page-003.png` col. 2 (`passive systems`, `A = 0`); value `"0"`, unit `pN / um^2`, `from_source: ref-23`, and quote exact.
- `passive_flat_growth_rate_bound` (p. 3): verified on `pages/page-003.png` col. 2 (`entropic cost prevents the aggregation of CMC`); value `"< 0"`, unit `1 / s`, `from_source: ref-23`, and quote exact.
- `protrusion_short_wavelength_growth_rate_sign` (p. 5): verified on `pages/page-005.png` col. 1 (`growth rate ω < 0`); value `"< 0"`, unit `1 / s`, `from_source: ref-34`, and quote exact.
- `protrusion_unstable_growth_rate_sign` (p. 5): verified on `pages/page-005.png` col. 1 (`initial perturbation grows (ω > 0)`); value `"> 0"`, unit `1 / s`, `from_source: ref-34`, and quote exact.
- `passive_cylinder_interaction_zero` (p. 6): verified on `pages/page-006.png` col. 1 (`J = 0 in equation (2.1)`); value `"0"`, unit `pN / um`, `from_source: ref-88`, and quote exact.
- `fig3c_time_i` (p. 6): verified on `pages/page-006.png` Fig. 3(c)(i) (`t = 0.002`); digits, unit, `from_source: ref-88`, and quote exact.
- `fig3c_time_ii` (p. 6): verified on `pages/page-006.png` Fig. 3(c)(ii) (`t = 0.1`); digits, unit, `from_source: ref-88`, and quote exact.
- `fig3c_time_iii` (p. 6): verified on `pages/page-006.png` Fig. 3(c)(iii) (`t = 1`); digits, unit, `from_source: ref-88`, and quote exact.
- `fig3c_time_iv` (p. 6): verified on `pages/page-006.png` Fig. 3(c)(iv) (`t = 2`); digits, unit, `from_source: ref-88`, and quote exact.
- `critical_tube_wavelength_ratio` (p. 7): verified on `pages/page-007.png` col. 1 (`wavelength that is twice the tube length [87]`); value `"2"`, unit `dimensionless`, `from_source: ref-87`, and quote exact.
- `concave_protrusive_growth_rate_sign` (p. 7): verified on `pages/page-007.png` col. 2 (`decaying sideways motion (damped waves)`); value `"< 0"`, unit `1 / s`, `from_source: ref-34`, and quote exact.
- `myosin_strength_zero_fig4c` (p. 8): verified on `pages/page-008.png` Fig. 4(c) caption (`The case of A* = 0`); value `"0"`, unit `pN / um^2`, `from_source: ref-116`, and quote exact.
- `fig4d_time_t1` (p. 8): verified on `pages/page-008.png` Fig. 4(d) caption (`t = 0, 1, 2, 3 s`); value `"1"`, unit `s`, `from_source: ref-116`, and quote exact.
- `fig4d_time_t2` (p. 8): verified on `pages/page-008.png` Fig. 4(d) caption (`t = 0, 1, 2, 3 s`); value `"2"`, unit `s`, `from_source: ref-116`, and quote exact.
- `fig4d_time_t3` (p. 8): verified on `pages/page-008.png` Fig. 4(d) caption (`t = 0, 1, 2, 3 s`); value `"3"`, unit `s`, `from_source: ref-116`, and quote exact.
- `fig4d_inset_time_t60` (p. 8): verified on `pages/page-008.png` Fig. 4(d) caption (`t = 60, 90, 120 s`); value `"60"`, unit `s`, `from_source: ref-116`, and quote exact.
- `fig4d_inset_time_t90` (p. 8): verified on `pages/page-008.png` Fig. 4(d) caption (`t = 60, 90, 120 s`); value `"90"`, unit `s`, `from_source: ref-116`, and quote exact.
- `fig4d_inset_time_t120` (p. 8): verified on `pages/page-008.png` Fig. 4(d) caption (`t = 60, 90, 120 s`); value `"120"`, unit `s`, `from_source: ref-116`, and quote exact.
- `myosin_wave_actin_lead_sign` (p. 8): verified on `pages/page-008.png` col. 1 (`actin concentration is at the leading front of the wave, while the myosin concentration is trailing`); value `"> 0"`, unit `um`, `from_source: ref-116`, and quote exact.
- `fig5d_wave_velocity_peleg2011` (p. 9): verified on `pages/page-009.png` Fig. 5(d) (`V = 0.67 μm/min`, correcting OCR `mm/min`); digits, unit `um / min`, `from_source: ref-41`, and quote exact (origin `reported` noted in `q_audit_qualitative_inequalities_no_verdict`).
- Unrecorded micrograph scale bars on `pages/page-004.png` (Fig. 2e `10 μm`, Fig. 2f `4.0 μm`): confirmed appropriately out of scope as stated in `source.yaml` notes.

## 2. Model & Equation Audit (`models/*` vs. page images & independent derivations)

- `free_energy_density_eq2_1` (p. 2): verified against Eq. (2.1) on `pages/page-002.png`; all 5 terms, signs, factors, natural log (consistent with Eq. 2.5), units, `run()`, and `symbolic()` exact.
- `bending_energy_density_eq2_1` (p. 2): verified against first term of Eq. (2.1) on `pages/page-002.png`; `run()` and `symbolic()` exact.
- `bending_energy_density_eq2_2` (p. 2): verified against Eq. (2.2) on `pages/page-002.png`; `(1/2)*(kappa*(1 - n) + kappa'*n)*(H - n*H_0)^2`, `run()`, and `symbolic()` exact.
- `membrane_height_dynamics_eq2_3` (p. 2): verified against Eq. (2.3) on `pages/page-002.png`; `xi*dh/dt = -delta F/delta h + A*(n - <n>)`, units, `run()`, and `symbolic()` exact.
- `cmc_total_current_eq2_4` (p. 3): verified against Eq. (2.4) on `pages/page-003.png` and by independent functional differentiation of Eq. (2.1); `J_tot = -(Lambda/n_s)*n*nabla(delta F/delta n) = J_diff + J_disp + J_agg + J_curv` exact.
- `cmc_diffusive_current_eq2_5` (p. 3): verified against Eq. (2.5) on `pages/page-003.png`; `J_diff = -D*nabla n` exact.
- `cmc_diffusion_einstein` (p. 3): verified against text below Eq. (2.6) on `pages/page-003.png`; `D = Lambda*k_B*T` exact.
- `cmc_dispersive_current_eq2_6` (p. 3): verified against Eq. (2.6) on `pages/page-003.png`; `J_disp = -(kappa*Lambda*H_0^2/n_s)*n*nabla n` exact.
- `cmc_aggregation_current_eq2_7` (p. 3): verified against Eq. (2.7) on `pages/page-003.png` and functional derivative of `-J*n^2 + (J/n_s)*(nabla n)^2` in Eq. (2.1); `(2*J*Lambda/n_s)*n*(nabla n + (1/n_s)*nabla(nabla^2 n))` exact (produces stabilizing `-q^4/n_s` term under `-nabla . J_agg`).
- `cmc_curvature_current_eq2_8` (p. 3): verified against Eq. (2.8) on `pages/page-003.png`; `(kappa*Lambda*H_0/n_s)*n*nabla H` exact.
- `cmc_steady_state_enrichment_eq2_9` (p. 3): verified against Eq. (2.9) on `pages/page-003.png` and steady state of `J_diff + J_curv = 0`; `prefactor: unknown` properly reflects printed `~`.
- `cmc_steady_state_with_dispersion` (p. 3): reconstructed Lambert-W steady state of `J_diff + J_disp + J_curv = 0` verified independently by direct integration; exact.
- `cmc_molecular_scales` (p. 3): checked against `a = 1/H_0 = 1/sqrt(n_s)` and `R = 1/H` on `pages/page-003.png`; `model.py` adds `abs()` in `run()` for negative convex `H_0, H` (logged as reading note `q_audit_molecular_scales_abs_convention`).
- `cmc_sorting_exponent_molecular` (p. 3): verified `chi = (kappa/(k_B*T))*(a/R)` and `enrichment = exp(chi)` against p. 3 col. 2; sign convention noted in `q_audit_molecular_scales_abs_convention`.
- `wavelength_from_wavenumber` (p. 4): verified against Fig. 2(a),(b) (`lambda = 2*pi/q`) on `pages/page-004.png`; exact.
- `flat_membrane_dispersion` (p. 3–5): independently re-derived 2x2 Jacobian of Eqs. (2.1), (2.3), (2.4) in Monge gauge `H = +nabla^2 h` (verified against Fig. 1a-c sign convention); all 4 entries `Gh, Gn, m12, m21` exact. (Note: `symbolic()` returns complex root in damped-wave regime `disc < 0`, matching `run()`'s real part in tested positive-feedback regimes.)
- `flat_membrane_dispersion_quasisteady` (p. 3–5): independently re-derived `xi -> 0` limit of `flat_membrane_dispersion`; exact algebraic cancellation of `kappa^2*H_0^2*n_0*q^2` confirmed.
- `flat_membrane_instability_thresholds` (p. 4–5): independently re-derived `A_c, sigma_c, q_c, lambda_c` from `det(M_2) = 0` at `J = 0`; all 4 expressions exact (assuming positive-feedback pairing `A*H_0 < 0`).
- `flat_membrane_fastest_mode` (p. 5): verified numerical maximization of `omega(q)` and `lambda_max = 2*pi/q_max`; exact (`q_max = 4.6675 um^-1`, `omega_max = 0.05964 s^-1`).
- `flat_membrane_passive_spinodal` (p. 3–4): verified `q -> 0` spinodal expressions `J_c(sigma > 0, q -> 0)` and `J_c(sigma = 0, q -> 0)`; finite-`q` caveat for small `sigma > 0` logged in `q_audit_passive_spinodal_finite_q`.
- `flat_membrane_protrusion_pde_1d` (p. 5): reconstructed IMEX-RK2 pseudospectral PDE solver verified term-by-term against Eqs. (2.1), (2.3), (2.4); linear/nonlinear splitting, dealiasing, positivity (`n_min > 0`), and `mode4_growth_rate = 0.059639` vs `omega_max = 0.059642` all verified across `t_end = 30, 50, 60 s`.
- `closed_contour_cell_shape_kabaso2011` (p. 5): verified `runnable: false` — closed-contour equations and parameters from [35] are absent from Gov (2018); linked to `q_closed_contour_model_unstated`.
- `tube_equilibrium_radius` (p. 6): verified `R_tube = sqrt(kappa/(2*sigma))` by minimizing bare Helfrich cylinder energy `2*pi*R*(kappa/(2*R^2) + sigma)`; exact.
- `cylindrical_ring_dispersion` (p. 6): **disagreement logged (`q_audit_cylinder_ctheta_geometry`)** — SymPy check proves `H_hh = (kappa/R_0^4)*(1 + (q*R_0)^2*(1 - H_0*R_0*n_0) + (q*R_0)^4)` omits the `1/sqrt(1 + h'^2)` slope factor in `c_theta = 1/(R*sqrt(1 + h'^2))`; under exact geometry, bare cylinder has `1 + (q*R_0)^4` and passive `A = 0, J = 0` circumferential-only Eq. (2.1) is strictly stable for all `q`.
- `cylindrical_ring_critical_mode` (p. 7): inherits the `cylindrical_ring_dispersion` geometry approximation (logged under `q_audit_cylinder_ctheta_geometry`); root-finding and maximization implementation verified.
- `cylindrical_critical_tube_length` (p. 7): verified `from_words: true` translation of p. 7 col. 1 (`L_c = pi/q_c = lambda_c/2`, ratio `2`); exact.
- `membrane_height_myosin_eq5_1` (p. 7): verified against Eq. (5.1) on `pages/page-007.png`; `xi*dh/dt = -delta F/delta h + A*(n - <n>) - A**(m - <m>)` exact.
- `myosin_kinetics_eq5_2` (p. 8): verified against Eq. (5.2) on `pages/page-008.png`; `dm/dt = -k_off*m + k_on*n` exact.
- `myosin_steady_state_eq5_2` (p. 8): verified `<m> = (k_on/k_off)*<n>` from Eq. (5.2) steady state; exact.
- `myosin_wave_dispersion` (p. 8): independently re-derived 3x3 Jacobian `M_3(q)` from Eqs. (2.1), (2.4), (5.1), (5.2); all entries exact. Confirmed `c_phase` is weakly dependent on `A` (~11% change over 4x change in `A` at `q = 2, xi = 10`), supporting p. 8 col. 1.
- `myosin_wave_dispersion_quasisteady` (p. 8): independently re-derived `xi -> 0, J = 0` quadratic `omega^2 + b*omega + c_0 = 0`, `c_ana = sqrt(kappa*Lambda*|H_0|*n_0*k_on*A*/(sigma*n_s))`, and `q_omega_approx = k_off/(2*c_ana)`; all exact (`c_ana` strictly independent of `A`).
- `myosin_acoustic_dispersion` (p. 8): verified against `omega''(q) proportional to q*c` on `pages/page-008.png` col. 1; `prefactor: unknown` exact.
- `myosin_wave_speed_scaling` (p. 8): verified against `c proportional to sqrt(A*)` on `pages/page-008.png` col. 1; `prefactor: unknown` exact.
- `myosin_wave_pde_1d_linear` (p. 8): linearized 3x3 Crank-Nicolson PDE matrix verified exact, **disagreement logged (`q_audit_wave_pde_speed_and_peak_tracking`)** on `wave_packet_speed = x_n / t_end` transient offset coincidence with `c_ana` and multi-crest `argmax` jump at `xi = 10`.
- `nonlocal_myosin_kernel_dispersion` (p. 8): verified `runnable: false` — non-local myosin kernel from [116] omitted in Gov (2018); linked to `q_nonlocal_myosin_kernel_unstated`.
- `two_curvature_wave_model_peleg2011` (p. 8–9): verified `runnable: false` — two-species PDEs and parameters from [41] omitted in Gov (2018); linked to `q_two_curvature_model_unstated`.
- `calcium_feedback_model_veksler2009` (p. 9): verified `runnable: false` — only feedback topology Fig. 6 from [123] printed; linked to `q_calcium_feedback_model_unstated`.

## 3. Unstated Inputs (`added.yaml`), Metadata, Relations, Refinements & Quotes Audit

- `added.yaml` (all 34 entries: `kBT_added`, `kappa_added`, `kappa_cmc_added`, `sigma_added`, `ns_added`, `H0_convex_added`, `H0_concave_added`, `H_local_convex_added`, `H_grad_added`, `n0_added`, `n_local_added`, `n_grad_added`, `n_grad3_added`, `Lambda_added`, `J_attr_added`, `xi_added`, `A_protrusive_added`, `A_contractile_added`, `A_star_wave_added`, `k_on_added`, `k_off_added`, `m_local_added`, `F_passive_added`, `q_unstable_probe_added`, `q_short_probe_added`, `q_cyl_probe_added`, `R_tube_cyl_added`, `sigma_cyl_added`, `L_protrusion_pde_added`, `n_pert_amp_pde_added`, `t_end_protrusion_pde_added`, `L_wave_pde_added`, `h0_dip_wave_pde_added`, `w0_dip_wave_pde_added`, `t_wave_pde_added`, `n0_cyl_passive_added`; pp. 2–8): every entry has a documented biophysical basis, valid units, and plausible alternatives (no exact identities given non-physical alternatives). Two behavioral sensitivities under alternatives were identified and logged: `sigma_added = 10.0 > sigma_c = 6.00 pN/um` stabilizes the flat protrusion probe (`q_audit_qualitative_inequalities_no_verdict`), and `xi_added = 10 pN*s/um^3` triggers multi-crest `argmax` switching in `myosin_wave_pde_1d_linear` (`q_audit_wave_pde_speed_and_peak_tracking`); `n0_cyl_passive_added` instability depends on the `c_theta` approximation (`q_audit_cylinder_ctheta_geometry`).
- `source.yaml` (pp. 1–13): title, author (`N. S. Gov`), year (`2018`), venue (`Phil. Trans. R. Soc. B 373: 20170115`), DOI (`10.1098/rstb.2017.0115`), `page_range: [1, 13]`, and printed page labels `2–13` (p. 1 has no printed folio) verified exact against `pages/page-001.png`–`pages/page-013.png`.
- `quantities.yaml` (pp. 2–9): all 51 quantities checked for pint dimensionality and physical distinctness; noted 1 unused quantity (`passive_dispersion_determinant_coeff`), retained per non-deletion rule.
- `references.yaml` (pp. 10–13): all 43 listed references cross-checked against `pages/page-010.png`–`pages/page-013.png`; every `from_source` and `described_from` id resolves cleanly with accurate citation text and DOIs.
- `relations.yaml` (all 14 relations, pp. 2–8): verified independently via `.venv/bin/boxgraph relations corpora/actin_bundling/records --only gov2018` — all 14 relations pass (`verified`).
- `refinements.yaml` (all 7 refinements, pp. 5, 8): checked dt and nx convergence on both maxima (`n_max`, `actin_myosin_lead`) and minima (`n_min`, `dn_min`) per `careful_v3.md`.
- Citation quotes (all `data.yaml`, `models/*/model.yaml`, `questions.yaml`, `relations.yaml`, pp. 2–9): automated single-column verification confirmed 0 cross-column interleaved quotes and 0 missing quotes across the entire record.

## 4. Fixes & Additions Made

- No genuine transcription errors were found in existing `data.yaml` numbers/signs/exponents/units or existing printed-equation `models/*/model.py` implementations.
- `J_diff_positive_grad_sign` & `cmp_J_diff_negative_sign` (p. 3): added missing printed sign prediction (`< 0` for `nabla n > 0` below Eq. 2.6) and reproduction comparison against `J_diff_eq2_5_val`.
- `J_disp_positive_grad_sign` & `cmp_J_disp_negative_sign` (p. 3): added missing printed sign prediction (`< 0` for `nabla n > 0` below Eq. 2.6) and reproduction comparison against `J_disp_eq2_6_val`.
- `concave_protrusive_oscillation_freq_sign` & `cmp_concave_protrusive_damped_wave_freq` (p. 7): added missing printed oscillatory frequency prediction (`> 0` for "damped waves" in Section 5) and reproduction comparison against `flat_disp_concave_protrusive_imag`.
- `myosin_wave_growth_rate_sign` & `cmp_myosin_wave_growth_sign` (p. 7–8): added missing printed wave-instability amplitude growth prediction (`omega' > 0`) and reproduction comparison against `myosin_disp_wave_real`.
- `myosin_wave_oscillation_freq_sign` & `cmp_myosin_wave_freq_sign` (p. 7–8): added missing printed wave-instability oscillatory frequency prediction (`omega'' > 0`) and reproduction comparison against `myosin_disp_wave_imag`.
- `fig4d_time_t0` (p. 8): added printed initial snapshot time `t = 0 s` from Fig. 4(d) caption (`t = 0, 1, 2, 3 s`), completing the reported snapshot time series.
- `run_flat_disp_at_qc`, `run_flat_disp_qs_at_qc`, & `cmp_flat_disp_qc_consistency` (p. 4–5): added chained runs evaluating full 2x2 and quasi-steady flat-membrane dispersion at `q = flat_qc_convex` and consistency comparison verifying identical neutral stability `omega(q_c) = omega_qs(q_c) = 0`.
- `cmp_myosin_acoustic_vs_qs_wave` (p. 8): added approximation comparison between `myosin_acoustic_omega_imag` (`q*c_ana`) and `myosin_disp_qs_wave_imag` (`omega''_qs(q)`) at `q = 4.0 um^-1`.

## 5. Audit Disagreements Logged (`kind: reading` in `questions.yaml`)

- `q_audit_cylinder_ctheta_geometry` (p. 6): disagreement with reconstructed `cylindrical_ring_dispersion` / `cylindrical_ring_critical_mode` / `q_cylinder_passive_condition_coverage` — omitting `1/sqrt(1 + h'^2)` in `c_theta` introduces a spurious `(q*R_0)^2*(1 - H_0*R_0*n_0)` term responsible for the passive `J = 0` instability (which vanishes identically under exact circumferential-only Eq. 2.1 geometry), and runs pass `H0_concave_added` where Gov p. 6 / Fig. 3(a),(c) label inner-surface circumferential ring CMC as convex.
- `q_audit_wave_pde_speed_and_peak_tracking` (p. 8): disagreement with reconstructed `myosin_wave_pde_1d_linear` / `cmp_pde_wave_speed_vs_c_ana` — `x_n / t_end` agrees with `c_ana` at `t_end = 10 s` only due to an initial Gaussian shoulder offset (`x_0 ~ 0.56 um`) masking the finite-`xi` propagation speed `0.115 um/s`, and global `argmax` jumps crests at `xi_added = 10`, yielding spurious negative `actin_myosin_lead = -1.26 um`.
- `q_audit_passive_spinodal_finite_q` (p. 4): refinement of `flat_membrane_passive_spinodal` / `q_passive_Tc_claim` — while `2*J*n_0 <= k_B*T*n_s` indeed guarantees passive stability for all `q > 0`, a tensioned membrane (`sigma > 0`) can undergo finite-`q` spinodal instability strictly below the `q -> 0` bound `J_c(sigma > 0, q -> 0)`.
- `q_audit_qualitative_inequalities_no_verdict` (p. 5, 9): reading note on verdict coverage (`< 0` / `> 0` sign comparisons receive `no verdict` so sign flips under alternatives like `sigma_added = 10.0` or `xi_added = 10` are not flagged automatically) and `origin: reported` tagging of `fig5d_wave_velocity_peleg2011`.
- `q_audit_molecular_scales_abs_convention` (p. 3): reading note on `abs()` added in `cmc_molecular_scales` (not printed on p. 3), decimal conversion `10^-3 -> 0.001` in `curvature_radius_ratio_range`, and sign dropping in `a/R = |H|/|H_0|`.

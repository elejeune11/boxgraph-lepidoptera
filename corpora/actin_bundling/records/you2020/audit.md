# Independent Audit Log: `you2020`

Audited record: `corpora/actin_bundling/records/you2020/`
Source PDF: `corpus/actin_bundling/2020_You.pdf` (You, Baskaran & Marchetti, *PNAS* 117, 19767–19772, 2020 + 13-page Supplemental Material and movie captions, PDF pp. 1–19).

All 19 pages of extracted text (`text/page-001.txt`–`page-019.txt`) and all 19 page images (`pages/page-001.png`–`page-019.png`) were inspected independently against `source.yaml`, `quantities.yaml`, `data.yaml`, `added.yaml`, `references.yaml`, all 21 models (`models/*/model.yaml` and `model.py`), `runs.yaml`, `comparisons.yaml`, `relations.yaml`, `refinements.yaml`, and `questions.yaml`.

---

## 1. Corrections to transcription / factual typos

- `source.yaml` (p. 14 / SI-8, `notes` line 76): Fixed numerical typo `v = 0.01939` -> `v = 0.01919` for the bi-supercritical 1D PDE speed at $\delta/\kappa = 4$ (computed `v_bisup_pde_dk4_pred = 0.019189`, refined to `0.019196` at $\Delta t = 1.0$, already recorded as `0.01919` in `questions.yaml` and `comparisons.yaml`).
- `data.yaml` (`rhoA0_0625`, p. 17 / SI-11): Fixed off-by-one row index in the `note` field from `"9th of 17 evenly spaced rows in [0.4, 0.8]"` to `"10th of 17 evenly spaced rows in [0.4, 0.8], 0-indexed row 9"` (rows step by `0.025` from `0.400` at row 1 to `0.600` at row 9 and `0.625` at row 10; the numeric value `"0.625"` was already correct).

---

## 2. Added printed results, runs, and comparisons

Previously uncompared predictions explicitly stated or plotted in the paper were added without deleting or altering any existing entry:

- `v_ratio_zero_printed`, `run_vratio_fs_dk1`, `cmp_v_ratio_fs_dk1` (p. 4 / 19770, Eq. 9 & Fig. 3B): Added exact sign-reversal threshold $v_A^s / v_B^s = 0$ at $\delta = \kappa$ separating blue ($< 0$) and green ($> 0$) branches (`reproduced`, `0.0` vs `"0"`).
- `im_lambda_h_real_zero_printed`, `cmp_figS2_im_dk0` (p. 9 / SI-3, Eq. S13 & Fig. S2): Added exact zero imaginary part $\mathrm{Im}[\lambda_\pm^h] = 0$ in the real-eigenvalue regime at $\delta = 0, \chi_A = -0.06$ (`reproduced`, `0.0` vs `"0"`).
- `cmp_ep_lambda1_plus_at_dc` and `cmp_ep_lambda1_minus_supercrit` (p. 12 / SI-6, Eq. S28a & Fig. S3a): Compared already-computed $\lambda_1^+(\delta = \delta_c) = 0$ and supercritical Goldstone branch $\lambda_1^-(\delta = 2\delta_c) = 0$ against `lambda1_goldstone_zero_printed` (`reproduced`).
- `u1_plus_3_figS3c_2dc_eye`, `cmp_figS3c_u1p3_2dc` (p. 12 / SI-6, Fig. S3c): Added eye-read third eigenvector component $u_1^+(3) \approx -0.035$ at $\delta/\delta_c = 2$ (completing Fig. S3c alongside panels a–b).
- `cmp_v_static_pde2d` (pp. 12–13 / SI-6–SI-7, Fig. S4g): Compared already-computed 2D static band velocity `v_pde2d_dk1_pred = 0` against `v_static_zero_printed` (`reproduced`).
- `cmp_figS4g_v_1mode_vs_2d` (pp. 12–13 / SI-6–SI-7, Fig. S4g): Added `approximation` comparison checking the SI Sec. IV claim that the 1D formula $v = \sqrt{\delta^2 - \delta_c^2}$ fits the 2D traveling band velocity at $\delta = 3\kappa$.
- `freq_2d_static_zero_printed`, `cmp_figS4h_freq_static_zero` (pp. 12–13 / SI-6–SI-7, Fig. S4h): Added exact zero oscillation frequency in the 2D static state at $\delta = \kappa < \delta_c$ (`reproduced`).
- `cmp_v_static_bisup_pde_dk1` (p. 14 / SI-8, Fig. S5c): Compared already-computed bi-supercritical 1D PDE static speed `v_bisup_pde_dk1_pred = 0` against `v_static_zero_printed` (`reproduced`).
- `cmp_fig3c_psit_ode` and `cmp_fig3c_ampratio_psit_ode` (p. 4 / 19770, Fig. 3C): Compared already-computed one-mode ODE final steady phase shift $\psi = 0.4636$ and amplitude ratio $\rho_A/\rho_B = 0.4472$ at $\delta/\kappa = 1.5$ against Fig. 3C.
- `run_fs_stab_chi006_dk2` (p. 10 / SI-4, Eq. S15): Added run of `static_fp_stability_eqS14_S16` at already-printed parameters $\chi_A = -0.06, \delta = 2\kappa$ (inside the $F_S$ existence domain, $R = (\rho_A^s)^2 = 0.005 > 0$) showing $\mathrm{Re}[\lambda_\pm^s] = 0$ (boundary of positive-real-part amplitude instability).
- `run_micro_mf_BA_printed` and `run_micro_mf_BA_pi12` (p. 16 / SI-10): Added BA cross-diffusivity runs alongside AB (`k_BA = 1, sigma_BA = 3`) confirming both $R_{AB} > 0$ and $R_{BA} > 0$ under the printed low-density mean-field formula.

---

## 3. Newly identified open questions and audit disagreements/nuances

- `q_eqS15_amplitude_eigenvalues_positive_real_part` (`possible_error`, p. 10 / SI-4, Eq. S15): SI Sec. II claims the coupled amplitude eigenvalues $\lambda_\pm^s$ *"always have a negative real part in the region of parameters where $F_S$ exists, hence are stable."* Exact algebra on Eq. (S14) gives $\det = 2\alpha_B R > 0$ and $\mathrm{tr} = -3R - \alpha_A - \alpha_B$ with $R = (\kappa^2-\delta^2)/\alpha_B - \alpha_A > 0$, so whenever $\alpha_A < -\alpha_B$ and $\kappa^2 + \alpha_B(\alpha_B - 2\alpha_A)/3 < \delta^2 < \kappa^2 - \alpha_A\alpha_B$ (inside the $F_S$ existence band between the black and dashed blue lines in Fig. 1A), **both** $\lambda_\pm^s$ have strictly positive real part (e.g., $\mathrm{Re}[\lambda_\pm^s] = 0$ at $\chi_A = -0.06, \delta = 2\kappa$; $\mathrm{Re}[\lambda_\pm^s] = +0.001675 > 0$ at $\chi_A = -0.05, \delta = 1.7\kappa$). The printed claim holds only for $\delta \le \delta_c$, not across the full existence region of $F_S$.
- `q_eqS19_Ds_modes_never_complex_conjugate` (`possible_error`, p. 11 / SI-5, Eqs. S19–S23): SI Sec. III states that substituting $\phi_A^0 = \phi_A^s$ into $D_0$ yields modes that *"change from real to complex conjugates at $\delta = \delta_c = \sqrt{\kappa^2 + \alpha_B^2}$."* In fact, $\Delta(\phi_A^s) = [(\delta_c^2 - \delta^2)/\alpha_B]^2 \ge 0$ identically for all real $\delta$: $\Delta$ touches zero tangentially at the exceptional point $\delta = \delta_c$ and never becomes negative, so the eigenvalues of $D_s$ (Eqs. S21–S23) are purely real on both sides of $\delta_c$.
- `q_micro_meanfield_same_sign_cross_diffusivity` (`unexplained_difference`, p. 16 / SI-10): Evaluating the printed low-density mean-field formula $R_{\mu\nu} = (k_r\sigma_r^4 + k_{\mu\nu}\sigma_{\mu\nu}^4)\frac{\pi}{6}$ (or $\frac{\pi}{12}$) at the paper's own traveling-state particle parameters ($k_r = 2000, \sigma_r = 1, k_{AB} = -2, k_{BA} = 1, \sigma_{AB} = \sigma_{BA} = 3$) yields $k_r\sigma_r^4 = 2000 \gg |k_{AB}|\sigma_{AB}^4 = 162$, giving **same-sign positive** cross-diffusivities $R_{AB} = 962.37 > 0$ and $R_{BA} = 1089.61 > 0$ ($|\delta_{AB}|/\kappa_{AB} = 0.081 \ll 1$), contrary to the opposite-sign requirement $\kappa_{AB}\kappa_{BA} < 0$ of main-text Eq. (4) because the naive $g(r) = 1$ integral overweights the excluded-volume core $r < \sigma_r$.
- `q_ap_domain_length_L64_vs_64p8_boundary` (`reading`, p. 17 / SI-11, `L_ap_added`): Audit nuance on the first reader's `L_ap_added = 64` (kept unchanged per audit rules). Across all $17 \times 13 = 221$ grid points of Fig. S7(a), $L = 64$ matches 220/221 linear onset classifications (including $(50, 0.625)$, $(80, 0.750)$, and $(80, 0.450)$), missing only $(v_0 = 67.5, \rho_A^0 = 0.475)$ where $L_{\min} = 64.64$; any $L \in (64.64, 65.00)$ matches all 221/221 points (upper-bounded by stable $(77.5, 0.45)$ at $L_{\min} = 65.00$ and $(50, 0.60)$ at $L_{\min} = 65.08$). Note also that the $t \in [0, 600]$ time scale of Fig. S7(c–d) scales as $L^2$, independently ruling out $L = 2\pi$ with rescaled $\kappa$.
- Disagreement/note on `cmp_bisup_dc_exact1m_vs_fitted` (p. 14 / SI-8): Although `dc_over_k_bisup_1m_pred = 1.6005` differs by 6.7% from the printed fit `"1.5"`, BoxGraph evaluates `cmp_bisup_dc_exact1m_vs_fitted` as `reproduced` solely because `kappa_0005` is printed to one significant digit (`"0.005"`), causing BoxGraph's input-rounding sensitivity band `[1.534, 1.667]` to overlap `[1.45, 1.55]`. No edit was made to `kappa_0005` because `"0.005"` is verbatim as printed; direct 1D PDE sweeps independently confirmed the bi-supercritical transition occurs between $\delta/\kappa = 1.58$ ($v = 0$) and $1.60$ ($v = 0.000632$), with $v = 0$ at $\delta/\kappa = 1.50$.
- Minor reader-choice observations (left unchanged per audit rules): `microscopic_ar_particles/model.yaml` (`runnable: false`, `stochastic: true`) omits `k_BB` (`= 1`) and `sigma_BA` (`= 3`) from its unused input port list and emits one expected validator warning because a non-runnable stochastic model cannot be refined; `run_vratio_at_dc` passes $\delta = \kappa$ amplitudes into its unused `v_ratio_general` port while evaluating `v_ratio_fs` at $\delta = \delta_c$.

---

## 4. Independent verification summary

- **All printed numbers in `data.yaml` (pp. 2–19)**: Checked every digit, sign, exponent, unit, page number, and quote against page images and text files; no numerical transcription or OCR exponent errors were found.
- **All 21 models (`models/*/model.yaml` + `model.py`)**: Independently re-derived and verified every printed equation — Eqs. (1)–(11) and SI Eqs. (S1)–(S40) plus unnumbered mean-field equations on p. 16 — symbol by symbol, sign by sign, and prefactor by prefactor. All 20 executable `model.py` implementations match the printed equations verbatim (with the $\pi/12$ correction properly placed in `variant_of: microscopic_ar_meanfield` model `microscopic_ar_meanfield_pi12`).
- **Independent re-verification of prior reader findings**:
  - Confirmed Fig. S1(a) (p. 8) rules out both the main-text subcritical-B model (`|phi_hat_B^1| = 0.2010`) and subcritical-B with $\phi_B^3$ restored (`|phi_hat_B^1| = 0.0865, |phi_hat_B^3| = 0.0126`), matching the bi-supercritical 1D PDE (`|phi_hat_A^1| = 0.1005, |phi_hat_B^1| = 0.1335, |phi_hat_A^3| = 0.00103, |phi_hat_B^3| = 0.00233`).
  - Confirmed SI Sec. V exact one-mode velocity $v = \sqrt{\delta^2 - \kappa^2 - (\kappa\alpha/\delta)^2}$, onset $\delta_c/\kappa = 1.6005$, and full PDE onset between $1.58$ and $1.60$.
  - Confirmed SI Sec. II numerical stability claim for $F_T$ (p. 10) across a broad parameter grid of Eq. (S11) (`max Re eig < 0` everywhere tested).
  - Confirmed SI Fig. S4(h) 2D oscillatory frequencies (`f = 0.00052, 0.00111, 0.00163, 0.00210` at $\delta/\kappa = 1.6, 2.0, 2.5, 3.0$) match the plotted curve while confirming the y-axis label `T` is a misprint for frequency.
  - Confirmed SI Sec. VI.A 2D polar integral gives $(\pi/12)\sigma^4$ (half of printed $\pi/6$) and white-noise variance misprint $\sqrt{2D}$ on p. 15.
  - Confirmed `added.yaml` contains all unstated inputs with clear bases and real verdict-changing alternatives, all 5 intra-record relations in `relations.yaml` pass `boxgraph relations`, and all 11 numerical refinements in `refinements.yaml` settle within tolerance.

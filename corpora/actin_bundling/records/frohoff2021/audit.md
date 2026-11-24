# Audit Log — `frohoff2021`

Independent audit of `corpora/actin_bundling/records/frohoff2021/` against `pages/page-001.png`–`pages/page-029.png` and `text/page-001.txt`–`text/page-029.txt` (Frohoff-Hülsmann, Wrembel & Thiele, *Phys. Rev. E* **103**, 042602, 2021; 29 pages including Appendices A–E).

## 1. Data, metadata, and references check (`data.yaml`, `source.yaml`, `references.yaml`)

- **`source.yaml`**: Verified title, authors (`Tobias Frohoff-Hülsmann`, `Jana Wrembel`, `Uwe Thiele`), venue (`Phys. Rev. E 103, 042602 (2021)`), DOI (`10.1103/PhysRevE.103.042602`), and page labels `1: "042602-1"` through `29: "042602-29"` against `pages/page-001.png`–`page-029.png`. No errors.
- **`data.yaml`**: Checked all 92 original data entries (65 reported parameters and 27 printed/derived predictions) against the rendered page images (`pages/page-001.png`–`page-026.png`) for digits, signs, exponents, approximate/inequality markers, units, `page`, `printed_page`, `origin`, `from_source`, and verbatim quotes in `text/page-NNN.txt`.
  - **0 transcription errors** found in `data.yaml`.
  - Minor note: `a_fig10b` quote (`"presents a space-time plot of a simulation at"`) wraps across two lines of the full-width Fig. 10 caption on `text/page-015.txt` (lines 49–50, no two-column interleaving; accepted by validator).
- **`references.yaml`**: Checked all 32 references against `pages/page-027.png`–`page-029.png`.
  - **Fixed 1 genuine transcription error**: `ref-53` (Brauns, Weyer, Halatek, Yoon & Frey, *Phys. Rev. Lett.* **126**, 104101) was transcribed as `(2020)` instead of `(2021)` as printed on `pages/page-028.png`. Corrected to `(2021)`.
- **Minor prose-note discrepancy (left untouched)**: In first-reader commentary (`source.yaml`, `comparisons.yaml`, `questions.yaml`, `refinements.yaml`), the lower three-phase hysteresis saddle-node fold in Fig. 5(b) is cited as `a = -1.0111`, whereas the actual converged spectral continuation (`run_bif_fig5b_hyst_dn`) produces `-1.010855` (refined to `-1.01085524` at `n_modes = 56`, i.e., `-1.0109`). Because this appears only in descriptive notes (and not in `data.yaml` or computed outputs), first-reader prose was not overwritten.

## 2. Model equations and numerical verification (`models/*`, `relations.yaml`, `refinements.yaml`)

Checked all 26 original models plus `coupling_path_eq23` (27 models total) against Eqs. (1)–(24), (A1)–(A8), (B1)–(B13), (C1)–(C3), (D1)–(D19), (E1)–(E3), and unnumbered inline formulas on `pages/page-003.png`–`page-027.png`:

- **0 equation, sign, or factor errors** found in any `model.py` implementation.
- **Dimensional consistency**: Verified `nondim_scales_eqA3_A6` (`tau`, `ell`, `Q`, `kappa`, `phi1_hat`, `phi2_hat`, `a`, `a_Delta`, `rho`, `alpha`) against Eqs. (A1)–(A8) on `pages/page-021.png`–`page-022.png` and dimensionalities in `quantities.yaml`.
- **Eq. (B12) misprint verification**: Independently inspected zoomed crop of `pages/page-022.png`. Confirmed Eq. (B12) row 1 prints `[q_-^o, q_+^o]` for $\kappa > 1$ (inverted since $1 - \kappa < 0 \implies q_+^o < q_-^o$), row 2 prints `[0, q_\mp^o]` for $\kappa \gtrless 1$ (correct), and row 3 prints `[q_-^o, q_+^o]` for $\kappa < 1$ (correct). Confirmed `complex_eigenvalue_band_eqB10_B12` faithfully implements the printed equation and `complex_eigenvalue_band_eqB10_B12_corrected` is a valid `variant_of` correction.
- **Independent re-verification of key first-reader findings**:
  - **Fig. 10(b) fastest mode (`a = -2.3`, p. 15)**: Independently evaluated Eq. (B1) for `n = 1..12`; confirmed `lam_+(n=8) = 30.9225 > lam_+(n=7) = 30.0189 > lam_+(n=9) = 27.1925 > lam_+(n=6) = 26.0590`, confirming `cmp_n_fastest_fig10b` (`not reproduced`) and `q_fig10b_fastest_mode_n6_vs_n8`.
  - **Fig. 6(a) / Fig. 7(b) `n = 3` stability at `a = 1.12` (p. 11)**: Independently built a full-domain `N = 144` spectral Jacobian (no Bloch reduction); top nonzero eigenvalues are `[-1.09286e-4, -1.09286e-4, -5.38684e-3, -5.38684e-3, ...]`, matching Bloch-Floquet `steady_state_bifurcation_eq7` identically and confirming continuum stabilization at `a = 1.121127` (`q_fig6a_n3_stability_a112`).
  - **Appendix C variational coexistence (`phi_1 = 2.5`, p. 23)**: Verified bulk free-energy Hessian eigenvalues are strictly positive for both coexisting phases A `(2.5, 1.2492)` (`[3.111, 17.321]`) and B `(2.3448, -1.2504)` (`[3.106, 15.079]`), confirming a thermodynamically stable binodal rather than a spinodal/unstable root.
  - **Refinements**: Verified all 10 entries in `refinements.yaml` settle cleanly (including extrema/minima convergence `ref_ss_fig7a_phi1_min` and `ref_pde_fig5d_ngrid_min` checking against Gibbs undershoot).
  - **Relations**: Verified all 8 original + 2 newly added relations in `relations.yaml` via `boxgraph relations` (all 10 verified symbolically or to `~1e-16`).

## 3. Additions made (integrating `coupling_path_eq23` and previously uncompared printed claims)

1. **Integrated `models/coupling_path_eq23` (Eq. 23, p. 14)**:
   - This model (added by an interrupted earlier audit pass and flagged `warning model coupling_path_eq23: never run`) was verified verbatim against Eq. (23) on `pages/page-014.png`: $\Delta(\Sigma) = \mp 2\rho\sqrt{\Sigma} + \Sigma$ for $\alpha \gtrless -\rho$, with single-column verbatim quote `"lows curves in Fig. 9 given by"`.
   - Added runs exercising both $\alpha > -\rho$ and $\alpha < -\rho$ branches in `runs.yaml`: `run_coupling_invariants_fig1c` + `run_coupling_path_fig1c` (`rho = 0.5, alpha = 0`), `run_coupling_path_fig8d` (`rho = 1.35, alpha = 1.65`), `run_coupling_invariants_fig13b` + `run_coupling_path_fig13b` (`rho = 1.3, alpha = -1.31`), and `run_coupling_invariants_fig9_rho03` + `run_coupling_path_fig9_rho03` (`rho = 0.3, alpha = -1.3`, exercising previously unused `rho_03`).
   - Added consistency comparison `cmp_Delta_path_fig1c` in `comparisons.yaml` comparing `Delta_path_fig1c_pred` against `Delta_fig1c = -0.25` (`reproduced`, exact).
   - Added stated relations `rel_eq23_same_as_invariants_pos` ($\alpha > -\rho$) and `rel_eq23_same_as_invariants_neg` ($\alpha < -\rho$) in `relations.yaml`, both verified by `boxgraph relations`.

2. **Fig. 6(a) / Fig. 7(b)–(c) `n = 3` stabilization claim (p. 11)**:
   - Previously `run_ss_fig6a_n3_a112` and `run_bif_fig6a_n3_stab` had no comparison in `comparisons.yaml`.
   - Added printed datum `a_fig6a_n3_stab_printed` (`"1.11 to 1.12"`, p. 11) and comparison `cmp_a_fig6a_n3_stab` against `a_fig6a_n3_stab_pred` (`1.12113` vs. `1.11 to 1.12` -> verdict `close`, quantitatively documenting the continuum shift above `1.12` in `q_fig6a_n3_stability_a112`).
   - Added `run_ss_fig6a_n3_a111` (at `a = 1.11, n = 3`, `warm_start: true`), datum `n_unstable_fig7c_printed` (`"0"`, p. 11), and comparison `cmp_n_unstable_fig7c` (`reproduced`, `0 == 0`).

3. **Fig. 12(c) standing-wave PDE simulation (p. 17)**:
   - `source.yaml` noted checking Fig. 12(a,c) and `a_fig12c = -0.59` was recorded in `data.yaml`, but no PDE run existed for Fig. 12(c).
   - Added `t_end_fig12c_added` (`"8000"`, alternatives `["6000", "10000"]`, matching Fig. 12(c) time axis `t x 10^3 = 0..8`) in `added.yaml`, `n_sw_fig12c_printed` (`"1"`) in `data.yaml`, `run_pde_fig12c` in `runs.yaml`, and `cmp_n_sw_fig12c` in `comparisons.yaml` (`reproduced`, `n_dominant_final = 1` with `norm_osc_std ≈ 0.105` across all seeds and alternatives).

4. **Additional runs exercising previously unused printed parameters**:
   - Added `run_subcrit_fig13b_n1` (`alpha = -1.31`, verifying `R_+ = +11.55 > 0` for the Fig. 13(b) `subplus` scenario of Eq. 24) and `run_subcrit_fig13c_n1` (`alpha = -1.315`, verifying `|alpha| > alpha_H = 1.31012` past the primary Hopf collision in Fig. 13(c)).
   - Added `run_ss_fig10b_n5` (`a = -2.3, n = 5`), confirming the paper's p. 15 claim that the steady `n = 5` state in Fig. 10(b) lies in the Hopf-unstable window (`lam_max = +0.00504`, `n_unstable = 4`).

## 4. Unstated inputs audit (`added.yaml`)

Audited all 42 original entries plus `t_end_fig12c_added` (43 total):
- Every entry has a documented `basis`, `stands_for`, and genuine physical or numerical `alternatives`.
- None attach non-physical alternatives to exact mathematical identities.
- All tested alternatives execute with zero failed runs and zero verdict changes.
- Minor note: `alpha_fig1ab_added` (`0.5`, alternatives `["1.0", "1.5"]`) varies $\Delta = \alpha^2 - \rho^2$ away from `0.25` at fixed `rho = 0`; this does not affect any comparison because the general Fig. 1 runs are uncompared.

## 5. Disagreements with first-reader choices (`reading` questions) and new open questions

Logged in `questions.yaml` without overwriting any first-reader questions or models:

- **`q_audit_hopf_endpoint_sign`** (`kind: reading`, p. 6): On p. 6 (below Eq. 21), the Hopf threshold endpoint on the CH curve is printed as $f_2''^{H_{\text{end}}} = \pm\sqrt{\Delta/Q}$ (with $-\sqrt{\Delta/Q}$ relevant on the stable boundary of Fig. 1b for $\kappa > 1/Q$ and $+\sqrt{\Delta/Q}$ in Fig. 1a for $\kappa < 1/Q$), whereas `hopf_threshold_general_eq14_15` prints and outputs only $+\sqrt{\Delta/Q}$.
- **`q_audit_drift_pitchfork_generalization`** (`kind: reading`, p. 17): On p. 17 (right column), the drift-pitchfork solvability condition is printed in its $Q = 1, \bar{\phi}_i = 0$ form $0 = \int (\phi_1^2 + \frac{\rho+\alpha}{\rho-\alpha}\phi_2^2)\,dx$, whereas `steady_state_bifurcation_eq7` and `bifurcation_threshold_finder_eq7` code the exact general form $\int [(\phi_1-\bar{\phi}_1)^2 + \frac{1}{Q}\frac{\rho+\alpha}{\rho-\alpha}(\phi_2-\bar{\phi}_2)^2]\,dx$ (identical on every run in this record).
- **`q_p17_decoupled_eigenvector_condition`** (`kind: possible_error`, p. 17): New paper-level finding — on p. 17 (left column, paragraph 1), for $\alpha \approx \rho$ where $\phi_2$ decouples from $\phi_1$, the paper states the eigenvector for $\lambda_2$ is $(-2k^2\rho/(\ell^2\lambda_1), 1)$. Direct linear algebra on Eq. (10) at $Q = 1, \alpha = \rho$ yields $(2k^2\rho/(\ell^2(\lambda_2 - \lambda_1)), 1)$, so the printed expression holds only at the bifurcation $\lambda_2 = 0$, a condition omitted in the prose.

## 6. Final validation status

- `.venv/bin/boxgraph validate corpora/actin_bundling/records/frohoff2021 --cache build/cache`: **SOUND (0 errors, 0 warnings)** across 27 models, 95 data entries, 43 added inputs, 64 runs, 38 comparisons (`21 reproduced`, `1 close`, `1 not reproduced`, `15 no verdict`), 10 refinements (all settled), and 0 alternative verdict changes.
- `.venv/bin/boxgraph relations corpora/actin_bundling/records --only frohoff2021`: **10/10 verified**.

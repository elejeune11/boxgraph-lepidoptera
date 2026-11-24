# Audit Report: `zwicker2015`

- **Source:** D. Zwicker, A. A. Hyman, and F. Jülicher, *"Suppression of Ostwald ripening in active emulsions"*, *Phys. Rev. E* **92**, 012317 (2015), DOI `10.1103/PhysRevE.92.012317`.
- **PDF audited:** `corpus/actin_bundling/2015_Zwicker.pdf` (PDF pages 1–13, printed labels `012317-1` to `012317-13`), inspected visually against every rendered page image `pages/page-001.png`–`pages/page-013.png` alongside `text/page-001.txt`–`text/page-013.txt`.
- **Scope:** Independent second-agent audit of `source.yaml`, `quantities.yaml`, `data.yaml`, `added.yaml`, all 63 models in `models/`, `runs.yaml`, `comparisons.yaml`, `refinements.yaml`, `relations.yaml`, `questions.yaml`, and `references.yaml`. (Figure digitization is left to the subsequent `digitize` job; no `figures.yaml` or `figures/` files were created.)
- **Baseline validator check:** `zwicker2015: SOUND (0 errors, 0 warnings)`.

---

## 1. Audit of `data.yaml`, `source.yaml`, and `references.yaml` Against Page Images

Every datum in `data.yaml` (41 original entries) was checked digit-by-digit, sign, exponent, unit, page number, printed page label, and quote against `pages/page-001.png`–`pages/page-013.png` and `text/page-001.txt`–`text/page-013.txt`:

- **Metadata (`source.yaml`):** Title, authors (`David Zwicker`, `Anthony A. Hyman`, `Frank Jülicher`), year (`2015`), venue (`Phys. Rev. E 92, 012317 (2015)`), DOI (`10.1103/PhysRevE.92.012317`), `page_range: [1, 13]`, and all 13 `012317-N` page labels match `pages/page-001.png`–`pages/page-013.png`.
- **References (`references.yaml`):** All 46 references (`ref-1` through `ref-46`) verified against `pages/page-013.png`; labels, authors, journals, volumes, pages, and years match.
- **Numeric transcription in `data.yaml`:** All 41 original numeric values, signs, powers of ten, units, and page numbers match the rendered page images:
  - Sec. III / Figs. 3–6 (pp. 4–7): `psi_binary = 1`, `kf_zero_passive = 0`, `kb_zero_passive = 0`, `phiB_bar_fig3a = 0.0909091` (printed as exact fraction `1/11`, converted to decimal because `boxgraph.precision` does not parse slash fractions), `phi0B_fig3a_target = 0.0909091`, `kf_over_k0_fig3a = 3.6 x 10^-5`, `kb_over_k0_fig3a = 3.6 x 10^-4`, `kb_over_kf_fig3b = 10`, `Vs_over_w3_fig3cd = 1.8 x 10^6`, `kb_over_kf_Bdroplet_bound = 3`, `kb_over_kf_symmetric = 1`, `kf_over_kb_Adroplet_bound = 3`, `phi0B_Bdroplet_bound = 0.25` (`1/4`), `phi0B_bicontinuous = ≈ 0.5` (`≈ 1/2`), `phi0B_Adroplet_bound = 0.75` (`3/4`), `Rbar_over_l_at_quarter = 1`, `kf_over_k0_fig4b = 1.4 x 10^-4`, `Vs_over_w3_fig4 = 5.8 x 10^4`, `fig4a_n_stable_states = 3`, `fig4a_n_small_unstable_states = 3`, `fig4a_n_large_sym_unstable = 1`, `fig4b_n_stable_states = 4`, `Vs_over_w3_fig5b_small = 1.5 x 10^5`, `Vs_over_w3_fig5b_large = 1.2 x 10^6`, `kb_over_k0_fig5c = 5 x 10^-3`, `kf_over_k0_fig5c = 5 x 10^-4`, `kf_over_kb_fig5d = 0.1`, `Vs_over_w3_fig6a = 2.3 x 10^5`.
  - Sec. IV / Figs. 7–10 (pp. 8–10): `kf_zero_autocat = 0`, `Q_zero_nocore = 0`, `psi_autocat = 0.1`, `k_over_k0_fig7 = 0.02`, `kb_over_k0_fig7 = 2 x 10^-4`, `phibar_fig7a = 0.02`, `Vs_over_w3_fig7 = 6.9 x 10^7`, `Q_over_w3k0_fig7 = 300`, `Q1_over_w3k0_fig10a = 1000`, `k_over_k0_fig9a = 100`, `kb_over_k0_fig9a = 0.01`, `Vs_over_w3_fig9a_small = 1.5 x 10^5`, `Vs_over_w3_fig9a_large = 8.8 x 10^5`, `fig8a_max_stable_droplets = 1`, `fig8b_n_stable_states = 1`, `fig8b_n_critical_radii = 0`, `fig10a_n_fixed_points = 1`.
  - Appendices E–F (p. 12): `pde_axial_points_appE = 128`, `pde_radial_points_appE = 32`, `homog_max_rate_sum_appF = 1`.
- **Quote fixes (genuine transcription/column-wrap fixes):** Automated single-line column verification across all citations detected 4 quotes spanning line breaks across adjacent columns in `pdftotext -layout`. All 4 were trimmed to exact single-column-line verbatim strings from `text/page-NNN.txt`:
  1. `data.yaml` (`fig8a_max_stable_droplets`, p. 8): `"at most a single droplet can be stable in this system."` → `"most a single droplet can be stable in this system."` (`page-008.txt` line 38).
  2. `models/binary_free_energy_sec3/model.yaml` (p. 4): `"The associated free energy density is f (φ B )"` → `"associated free energy density is f (φ B )"` (`page-004.txt` line 74).
  3. `questions.yaml` (`q_binary_free_energy_sec3_missing_psi2`, p. 4): same single-line quote fix (`page-004.txt` line 74).
  4. `models/unnormalize_binary_params/model.yaml` (p. 5): `"quantities are normalized to the length scale w ="` → `"quantities are normalized to the length scale"` (`page-005.txt` line 38).
  After these fixes, 100% of quotes across `data.yaml`, `models/*/model.yaml`, `questions.yaml`, and `relations.yaml` appear verbatim on a single line of their cited `text/page-NNN.txt`.

---

## 2. Audit of All 63 Models Against Printed Equations

Every model's `model.yaml`, `model.py` (`run()` and `symbolic()`), symbols, signs, factors, and units were audited against the equations printed on `pages/page-002.png`–`pages/page-013.png`, and every analytical derivation in Appendices A–G and Secs. II–IV was independently re-derived in SymPy:

- **Sec. II continuum & interface equations (pp. 2–4, Eqs. 1–12):**
  - `ternary_free_energy_eq1`, `flat_interface_profile`, `interface_width_w`, `surface_tension_gamma`, `laplace_coeff_beta`, `capillary_length_gb`, `capillary_length_from_micro`, `diffusivity_from_mobility`, `reaction_diffusion_length_l`, `reaction_flux_eq4`, `surface_flux_eq8_11`, `surface_flux_small_R`, `coarse_grained_droplet_eq9_12` all match the printed equations exactly.
  - Note: $D_A = m_A a$ (p. 3 below Eq. 7) is not modeled separately because $D_A$ never enters the coarse-grained droplet dynamics under the homogeneous-$\phi^A$ approximation ($D_A \approx D_B$).
  - Note on `gibbs_thomson_eq6`: implements Eq. (A2) ($\phi_-^B = \psi + \gamma\beta/R$, $\phi_+^B = \gamma\beta/R$) rather than the truncated Eq. (6) ($\phi_-^B \approx \psi$); logged as a minor reading note in `q_audit_reading_notes`.
- **Sec. III binary first-order emulsion equations (pp. 4–7, Eqs. 13–21, Fig. 6b, Sec. V scalings):**
  - `binary_free_energy_sec3` (printed without $1/\psi^2$), `binary_droplet_growth_eq13`, `binary_background_ss_eq15`, `binary_critical_radius_approx`, `binary_single_radius_eq16`, `binary_critical_rate_eq17`, `binary_critical_kb_fig6b`, `binary_single_droplet_ss_solver`, `binary_coarsening_rate_eq18`, `passive_coarsening_rate_eq18`, `binary_stability_radius_eq19`, `binary_dense_volume_eq20`, `binary_stability_number_eq21`, `binary_multidroplet_stability_solver`, `system_radius_Rs`, `interface_diffusion_rate_k0`, `ratio_Rbar_over_l`, `Rstab_kb_scaling`, `Rbar_kb_scaling`, `lifshitz_slyozov_scaling_printed`, `lifshitz_slyozov_scaling_corrected`, and processing helpers (`rate_from_ratio`, `forward_rate_from_ratio`, `unnormalize_binary_params`) all match the printed equations accurately.
  - Independently confirmed that at the saddle-node bifurcation $k_f = k_f^c$ (Eq. 17), the exact single-droplet radius $R_* = 3\gamma\beta/(2\phi_0^B)$ equals $R_{\rm stab}$ (Eq. 19) identically, proving the paper's p. 6 claim that $\bar{R} > R_{\rm stab}$ whenever $k_f < k_f^c$ in the low-density regime.
  - Independently confirmed the Lifshitz-Slyozov dimensional misprint on p. 7: printed $\langle R \rangle \propto (t D_B / w^2)^{1/3}$ is dimensionless on the RHS (missing factor $w$).
- **Sec. IV ternary autocatalytic equations (pp. 8–10, Eqs. 22–25):**
  - `autocatalytic_growth_eq22`, `effective_rate_ks`, `autocatalytic_volume_eq23`, `autocatalytic_critical_radius`, `autocatalytic_nucleation_Qc`, `ternary_coarsening_rate_eq24`, `ternary_coarsening_rate_eq24_corrected`, `autocatalytic_stability_Qstab_eq25`, `unequal_cores_volume_ratio`, `autocatalytic_emulsion_solver`, and `unnormalize_ternary_params` match the printed equations.
  - Independently confirmed the Eq. (24) misprint: substituting Eqs. (G2)–(G3) into Eq. (G1) yields $Q/(\bar{V} k_s)$ identically, whereas printed Eq. (24) omits subscript $s$ ($k$ instead of $k_s = k_f + k\psi$).
  - Independently confirmed `q_sec4_length_scale_w_inconsistency`: Figs. 7(a), 8(a), 8(b), 9(b), 10(a), and 10(b) all use capillary normalization $\gamma\beta = w/6$ (e.g., Fig. 7a $R_{\rm crit}/w = 16.7$, peak $\dot{R}/(wk_0) = 1.5\times 10^{-3}$, Fig. 9b red line $Q_{\rm stab}/(w^3 k_0) = 209.4$, Fig. 10a $(R_1, R_2)/w = (109.9, 63.4)$), whereas Fig. 9(a) uses the continuum relation $\gamma\beta = w\psi/6 = w/60$ ($\lambda/k_0 = 7.0\times 10^{-5}$ and $2.1\times 10^{-4}$ at $\bar{V}/w^3 = 10^4$).
  - Independently confirmed `q_Qc_vs_Qstab_inequality_tension`: Eq. (25) at $k_f = 0, \phi_0^B = 0$ gives $Q_{\rm stab} = 4\pi D_B \gamma\beta k / k_b$, identically equal to $Q_c$ as defined on p. 8 with $\phi_0^A = k_b/k$; the strict inequality $Q_{\rm stab} > Q_c$ for $k > k_b/\bar{\phi}$ (p. 9 col 2) requires evaluating $Q_c$ at the droplet-free trivial state $\phi_0^A = \bar{\phi}$.
  - Also verified Eq. (25) vs Appendix G: solving corrected Eq. (24) for $\lambda = 0$ gives $Q_{\rm stab}^{(24)} = 4\pi D_B k_s R(3\phi_+^B - 2\phi_0^B) / [3\psi(k_b + J/(\psi\bar{V}))]$, which reduces to Eq. (25) after additionally neglecting $J/(\psi\bar{V}) \ll k_b$; remarkably, solving exact Eq. (G1) without the $Q \ll k_s\bar{V}$ expansion gives $Q_{\rm stab}^{\rm exact} = 4\pi D_B k_s R^3(3\phi_+^B - 2\phi_0^B) / [3(R^2 k_b\psi - D_B\phi_0^B)]$, whose $\phi_0^B = 0$ limit matches Eq. (25) *exactly* with zero approximation error.
- **Appendices A–G (pp. 11–13):**
  - **Appendix A** (`coexistence_exact_eqA1`): SymPy first-order expansion confirms $\phi_-^B \approx \psi + \beta\gamma H$ and $\phi_+^B \approx \beta\gamma H$ with $\beta = 2/(b\psi)$.
  - **Appendix B** (`autocatalytic_volume_quad_appB`, `autocatalytic_volume_eqB1`, `autocatalytic_volume_eqB1_corrected`): SymPy confirms the exact quadratic $k\psi^2 V^2 - V_s\psi(k\bar{\phi} - k_b - Q/V_s)V - V_s(Q\bar{\phi} - 4\pi D_B\gamma\beta) = 0$ and that printed Eq. (B1) differs from its exact linear expansion by $-4\pi D_B\gamma\beta/(k_b\psi)$ (an extra factor $k\bar{\phi}/k_b$ on the surface-tension term).
  - **Appendix C** (`outer_profile_eqC1`): SymPy confirms Helmholtz residual = 0 and boundary flux reproduces Eq. (11).
  - **Appendix D** (`binary_jacobian_eqD3_D4`, `binary_jacobian_eqD3_D4_corrected`): SymPy confirms the spurious extra $N$ in printed Eq. (D3b) ($J_{0i} = -4\pi D_B N \bar{\phi}_0^B V_s^{-1}$ vs exact $-4\pi D_B \bar{\phi}_0^B V_s^{-1}$). **However, we disagree with the closing remark of `q_eqD3b_spurious_N_factor`** (see Section 5 below).
  - **Appendix E** (`binary_axisymmetric_pde_eq3b`): Nondimensionalization of Eq. (3b) matches exact Cahn-Hilliard dynamics ($\partial_{\tilde{t}}\phi = \tilde{\nabla}^2[\phi - 3\phi^2 + 2\phi^3 - \frac{1}{4}\tilde{\nabla}^2\phi] + \tilde{k}_f(1-\phi) - \tilde{k}_b\phi$). Numerical caveats logged in `q_audit_axi_pde_numerics` (see Section 5).
  - **Appendix F** (`homogeneous_dispersion_eqF1`, `_corrected`, `homogeneous_instability_eqF2`, `_corrected`, `binary_cahn_hilliard_1d_pde`): SymPy confirms exact linearization of Eq. (3b) has half the printed $q^2$ coefficient in Eq. (F1) (while matching the $q^4$ term $-D_B w^2 q^4 / 4$), reducing the instability cutoff in Eq. (F2) by a factor of 4 from $(k_f + k_b)/k_0 < 1$ to $(k_f + k_b)/k_0 < 1/4$, independently verified by `binary_cahn_hilliard_1d_pde` (`0.250000`).
  - **Appendix G** (`ternary_coarsening_rate_eqG1`): SymPy confirms Eqs. (G1)–(G3) yield `ternary_coarsening_rate_eq24_corrected` identically.

---

## 3. Genuine Wiring Errors Fixed & New Claims/Checks Added

### 3.1 Genuine Input Cross-Wiring Bug Fixed in `runs.yaml`
Inspecting all 72 run outputs revealed that 6 Sec. IV (ternary autocatalytic) runs in `runs.yaml` had been cross-wired to Sec. III (binary Fig. 3a) outputs (`phi0B_fig3a_pred = 1/11 = 0.09091` and `J_smallR_fig3a = -9.33e-18 m^3/s`), producing unphysical negative outputs:
- `run_coarse_grained_eq9_12_fig7` (which also combined single-droplet `Rbar_eq23_fig7` with `N_pair = 2`): produced negative `phi_0_A_cons = -0.0727` (since $\phi_0^B = 0.0909 > \bar{\phi} = 0.02$).
- `run_autocat_growth_eq22_fig7`: inherited negative `phi0A_eq12_fig7 = -0.0727` and binary `phi0B_fig3a_pred = 0.0909`.
- `run_Qstab_eq25_fig9b`: produced unphysical negative `Q_stab = -8.785e-15 m^3/s` (because $3\phi_+^B - 2\phi_0^B = 0.0042 - 0.1818 < 0$).
- `run_ternary_rate_eqG1_fig7`, `run_ternary_rate_eq24_printed_fig7`, `run_ternary_rate_eq24_corrected_fig7`: used negative binary Fig. 3a influx `J_smallR_fig3a` instead of positive Fig. 7 autocatalytic efflux.

**Fix made:**
- Added `phi0B_probe_autocat = 0.001` (alternatives `["0.0", "0.0013"]`, satisfying $0 \le \phi_0^B < \phi_+^B(\bar{R}) \approx 1.41\times 10^{-3}$) in `added.yaml`.
- Added `run_surface_flux_autocat_fig7` computing physical Fig. 7 efflux `J_smallR_autocat_fig7 = +6.103e-19 m^3/s`.
- Rewired the 6 Sec. IV runs to `phi0B_probe_autocat`, `J_smallR_autocat_fig7`, and `R1_pair_fig8b` (for `run_coarse_grained_eq9_12_fig7` with `N_pair = 2`). All outputs are now physically consistent (`phi_0_A_cons = +0.01006`, `Q_stab = +1.105e-16 m^3/s`, `lam_G1_fig7 = lam_24_corr_fig7 = -2.766e-6 1/s`), and printed Eq. (24) now predicts the *wrong sign* (`lam_24_printed_fig7 = +1.165e-6 1/s`, unstable at $Q = 300\,w^3 k_0$ where Fig. 8b is stable).

### 3.2 Additional Printed Claims, Refinements, and Relations Added
1. **Passive max stable droplets (p. 6 col 1, Fig. 4a):** Added datum `fig4a_max_stable_droplets = "1"`, exposed `N_max: fig4a_max_stable_pred` from `run_fig4a_multidroplet`, and added `cmp_fig4a_max_stable_droplets` (verdict: **reproduced**).
2. **Unequal catalytic cores with one core below $Q_{\rm stab}$ (p. 10 col 1, Fig. 10b):** Added datum `fig10b_one_below_Qstab_max_droplets = "1"`, added input `Q1_below_Qstab_fig10b = "100"` (alternatives `["50", "150"]`), added runs `run_params_fig10b_core100` and `run_fig10b_one_below_Qstab`, and added `cmp_fig10b_one_below_Qstab` (verdict: **reproduced** at $Q_1 = 50, 100$, but **depends on added input**: alternative $Q_1 = 150 < Q_{\rm stab} = 209.4$ yields `2` stable droplets inside the bistable window — see Section 5.2).
3. **Active axisymmetric PDE positivity refinement:** Added `ref_pde_axi_active_dt_phi_min` in `refinements.yaml` (using already-cached `dt` refinement series of `run_pde_axi_fig5c_active`).
4. **Intra-record relations added to `relations.yaml`:**
   - `rel_binary_free_energy_sec3_eq1_psi1` (`stated`, `same_as` at $\psi = 1, \phi^A = 0$): **verified** symbolically.
   - `rel_binary_free_energy_sec3_vs_eq1_general_psi` (`proposed`, `same_as` for $\psi \in [0.1, 0.8]$): **refuted** (factor $\psi^2$ discrepancy from missing $1/\psi^2$ in Sec. III).
   - `rel_eqD3d_corrected_to_eq18` (`stated`, `same_as` verifying `e_tilde` of `binary_jacobian_eqD3_D4_corrected` matches Eq. 18 identically): **verified** symbolically.
   - `rel_eq24_small_sys_to_passive_eq18` (`stated`, `same_as` verifying Eq. 24 at $Q = 0, J = 0, \psi = 1$ reduces to passive Eq. 18): **verified** symbolically.

---

## 4. Audit of `added.yaml`

All 26 original entries (plus 2 added during audit) have clear bases and genuine physical/numerical alternatives:
- Dimensional reference scales (`w_ref`, `DB_ref`, `b_ref`, `a_ref`, `kappa_ref`, `mB_ref`) span realistic orders of magnitude and leave all dimensionless comparisons invariant (`verdict_changes: []`).
- Reconstructed PDE solver inputs (`pde_1d_*`, `axi_*`) are properly exposed per `careful_v3` #5 rather than hardcoded.
- **Minor observation (`N_pair`):** `N_pair = 2` is placed in `added.yaml` with `alternatives: ["3", "5"]`, even though two droplets ($N = 2$) is explicitly stated in the captions of Figs. 4, 5, 8, and 10. Per audit rules, `N_pair` was not deleted; this classification note is recorded in `q_audit_reading_notes`.
- **Minor observation (`r_probe_outer` vs `R_probe_small`):** Under the alternative `R_probe_small = 5.0e-5`, `r_probe_outer = 2.0e-5 < R_probe_small` evaluates Eq. (C1) at $r < R$ without crashing (`failed_runs: []`, no comparisons affected).
- **Active sensitivity (`Q1_below_Qstab_fig10b`):** Alternative `"150"` flips `cmp_fig10b_one_below_Qstab` from `reproduced` (`1.0`) to `close` (`2.0`, within single-digit integer `1 ± 0.5` close band `[-0.5, 2.5]`), directly exposing the bistability finding in Section 5.2 with zero failed runs across all 56 alternative evaluations.

---

## 5. Open Questions & Audit Disagreements Recorded in `questions.yaml`

All 12 original questions in `questions.yaml` were reviewed and preserved without deletion. Four new questions were added:

1. **`q_audit_appD_exact_eigenvalues_always_real` (`kind: reading` — Disagreement with `q_eqD3b_spurious_N_factor`):**
   - The closing remark of `q_eqD3b_spurious_N_factor` claims that the paper's statement $e_- < e_+ < \tilde{e}$ (p. 12) ignores regimes where $(J_{00} - J_{ii})^2 + 4 N J_{0i} J_{i0} < 0$.
   - **Audit proof:** At every steady state of Eq. (13), $B \equiv D_B\bar{\phi}_0^B/\bar{R}^2 = D_B\gamma\beta/\bar{R}^3 + k_b/3$. Setting $A \equiv 4\pi D_B N \bar{R}/V_s > 0$ gives $J_{ii} = B - k_b$, $J_{00} = -(k_f + k_b) - A$, $J_{00} - J_{ii} = -(k_f + A + B) < 0$, and (for exact $J_{0i}$ without the misprinted $N$) $4 N J_{0i} J_{i0} = -4 A B$. Hence $\Delta_{\rm exact} = (k_f + A + B)^2 - 4 A B = (A - B)^2 + k_f^2 + 2 k_f(A + B) > 0$ **strictly for all physical steady states** (verified across 50,000 random steady states). Thus the paper's claim $e_- < e_+ < \tilde{e}$ (all real) is exact for Eqs. (13)–(14); complex conjugate pairs occur *only* because of the misprinted extra $N$ in Eq. (D3b) (which drives $\Delta_{\rm printed} < 0$ at Fig. 4b and across most $N \ge 2$ steady states) or off the steady-state manifold.

2. **`q_audit_sec4_autocat_bg_bistability` (`kind: overstated_claim` — New Finding on Sec. IV.B, Figs. 9b & 10b):**
   - Solving Eq. (22) coupled to full Eq. (10) (including $+k\phi_0^A\phi_0^B$) and Eq. (12) matches Figs. 7(a), 8(a), 8(b), 10(a), and 10(b) steady volumes ($V_2(Q_1=0) = 6.27\times 10^6 w^3$, Fig. 7a finite peak $\dot{R}/(wk_0) = 1.43\times 10^{-3}$ at $R/w = 32.7$).
   - However, because droplet steady state forces $k\phi_0^A \approx k_b$, the net background decay rate $k_b - k\phi_0^A$ in Eq. (10) cancels by ~60×, keeping $\phi_0^B \approx 0.93\,\phi_+^B$ instead of $\phi_0^B \approx 0$ even though $V_s/(N l^3) \approx 98 \gg 1$ with $l = (D_B/k_b)^{1/2}$ (in tension with p. 9 col 1).
   - As a result, linear Ostwald stability ($\lambda < 0$) of two droplets begins far below Eq. (25)'s $\phi_0^B = 0$ threshold ($Q_{\rm stab} = 209.4\,w^3 k_0$): at $Q \approx 79.0\,w^3 k_0$ in Fig. 9(b) ($\bar{\phi} = 0.02$) and at a saddle-node fold $Q_1 \approx 133.5\,w^3 k_0$ in Fig. 10(b) ($Q_2 = 300\,w^3 k_0$). Thus p. 10's claim that only a single droplet is stable whenever $Q_1 < Q_{\rm stab} < Q_2$ is overstated (bistable for $133.5 < Q_1/(w^3 k_0) < 213.9$). The sharp transitions plotted near $Q \approx 212\text{--}214\,w^3 k_0$ in Figs. 9(b) and 10(b) instead match the 2nd-droplet spontaneous nucleation threshold $Q\phi_0^A = 4\pi D_B\gamma\beta$ ($212.6$ in Fig. 9b; $213.9$ in Fig. 10b).

3. **`q_audit_axi_pde_numerics` (`kind: numerics` — Reconstructed 3D Axisymmetric PDE Caveats):**
   - `ref_pde_axi_active_dt_lam` (`run_pde_axi_fig5c_active`) exhibits a 4.6% shift from `dt = 1.0` (`-0.002383`) to `dt = 0.125` (`-0.002498`), exceeding its 2% tolerance (`belongs_to: numerics`); `dt <= 0.25` is required for <1% time-step convergence.
   - `run_pde_axi_fig5a_passive` initializes droplets with mean volume fraction `0.09327` (2.6% above the Fig. 5a target `1/11 = 0.09091`) because the code only adjusts upward when `curr < phi_bar`.
   - `phi_min` is initialized to `min(phi_0) = 0.0` and sampled every `nsteps // 30` steps, so `ref_pde_axi_passive_dt_phi_min` and `ref_pde_axi_active_dt_phi_min` report `0.0` identically (missing a negligible step-1 transient `-4.5e-6` in passive, while late-time `min(phi) = 0.0127 > 0` confirms positivity).

4. **`q_audit_reading_notes` (`kind: reading` — Minor Reading Notes):**
   - Records `N_pair` in `added.yaml` vs stated $N = 2$; `gibbs_thomson_eq6` using Eq. (A2) $\phi_-^B = \psi + \gamma\beta/R$ vs Eq. (6) $\phi_-^B \approx \psi$; `homog_max_rate_sum_appF` `"1"` scoring `0.25` as `close` due to single-digit precision width `±0.5` despite an exact factor-of-4 discrepancy; and `autocatalytic_emulsion_solver` raising `ValueError` at `Q_1 = 0, Q_2 > 0` (to flag for downstream Fig. 10b digitization).

---

## 6. Final Validation Summary

- Command: `.venv/bin/boxgraph validate corpora/actin_bundling/records/zwicker2015 --cache build/cache`
- Result: **`zwicker2015: SOUND (0 errors, 0 warnings)`**
- Record totals after audit: 43 data entries, 28 added inputs, 63 models, 75 runs, 18 comparisons (15 `reproduced`, 2 `close`, 1 `no verdict`; 1 verdict-changing alternative `Q1_below_Qstab_fig10b = 150`; 0 failed runs), 21 refinements, 31 relations (all 4 newly added relations verified/refuted as expected), 16 open questions, 46 references.
- Zero existing data, models, or questions deleted; strict writable path compliance (`corpora/actin_bundling/records/zwicker2015/**` only; no scratch files or `figures.yaml`/`figures/`).


# Figure Digitization Log — Zwicker, Hyman & Jülicher (2015)

Source: `corpus/actin_bundling/2015_Zwicker.pdf` (*Phys. Rev. E* **92**, 012317, PDF pages 1–13, Figs. 1–10).

## Overview

- **Total figures in paper**: 10 (Figs. 1–10 across pages 2–10; pages 1 and 11–13 contain text, equations, appendices, and references only).
- **Qualitative / non-quantitative panels left out**: 6 panels (Figs. 1a, 1b, 2a, 2b schematics; Fig. 5a and 5c top grayscale PDE density snapshots).
- **Quantitative panels digitized**: **22 panels** across Figs. 3–10, yielding **58 figure readings** (`figures.yaml`, `figures/*.py`) and **58 data entries** (`data.yaml`: 45 `predicted`, 11 `run`, 2 `fitted`).
- **Resolution & calibration**: Every quantitative panel was rendered at **1200 dpi** and calibrated automatically from sub-pixel black box spines, inward major/minor ticks, and (where present) horizontal zero lines via `figures/_calib.py`. Maximum axis calibration residual across all 22 panels is **0.428 px** (0.000–0.386 px for 21 of 22 panels).
- **Visual verification**: All 58 check overlays were rendered via `.venv/bin/boxgraph figure check corpora/actin_bundling/records/zwicker2015 <reading_id> --out corpora/actin_bundling/records/zwicker2015/figwork/<reading_id>.png` and visually inspected to confirm magenta crosshairs sit squarely on the plotted curves, boundaries, and markers.
- **Model testing**: 17 scale-invariant normalized wrapper models (preserving all 63 existing models unchanged), 61 swept runs (`runs.yaml`), and **87 figure comparisons** (`comparisons.yaml`: 54 `reproduced`, 8 `close`, 25 `not reproduced`) were added and validated across all printed vs. corrected `variant_of` models and length-scale conventions (`.venv/bin/boxgraph validate corpora/actin_bundling/records/zwicker2015 --cache build/cache` → `SOUND (0 errors, 0 warnings)`).

---

## Complete Figure & Panel Inventory (Figs. 1–10)

### Qualitative Figures & Panels Left Out

| Figure / Panel | Page | Description | Status | Reason Left Out |
|---|---|---|---|---|
| **Fig. 1(a)** | 2 | Schematic of full continuum fields $\phi^A, \phi^B$ with diffuse interfaces | Left out | Qualitative schematic without quantitative axes or data |
| **Fig. 1(b)** | 2 | Schematic of coarse-grained droplet radii $R_i$ and background fractions $\phi_0^A, \phi_0^B$ | Left out | Qualitative schematic without quantitative axes or data |
| **Fig. 2(a)** | 3 | Reaction network diagram for first-order $A \rightleftharpoons B$ kinetics | Left out | Qualitative reaction schematic |
| **Fig. 2(b)** | 3 | Reaction network diagram for autocatalytic $A + B \to 2B$, $B \to A$, and active core $Q_i$ | Left out | Qualitative reaction schematic |
| **Fig. 5(a) top strip** | 7 | Grayscale axisymmetric PDE $\phi^B(r,z)$ snapshots of two coarsening passive droplets | Left out | Qualitative density snapshots; droplet volumes $V_1(t), V_2(t)$ and $\Delta V(t)$ are quantitatively digitized from the lower panels |
| **Fig. 5(c) top strip** | 7 | Grayscale axisymmetric PDE $\phi^B(r,z)$ snapshots of two equalizing active droplets | Left out | Qualitative density snapshots; droplet volumes $V_1(t), V_2(t)$ and $\Delta V(t)$ are quantitatively digitized from the lower panels |

---

### Quantitative Panels Digitized (22 Panels, 58 Readings, All at 1200 dpi)

| Panel | Page | Crop `[x0, y0, x1, y1]` | Axis Scales & Calibration Misfit (`x_res`, `y_res`) | Reading IDs (`n` pts) | Key Digitized Values & Notes |
|---|---|---|---|---|---|
| **Fig. 3(a)** | 5 | `[800, 1200, 2900, 2620]` | lin / lin; `0.200 px`, `0.350 px` | `fig3a_passive` (12), `fig3a_active` (12) | Passive peak `0.01241` at `R/w = 3.67`; active crosses 0 at `R_crit/w ≈ 1.84`, `R_bar/w ≈ 26.5` |
| **Fig. 3(b)** | 5 | `[2920, 1200, 4960, 2620]` | log10 / log10; `0.337 px`, `0.266 px` | `fig3b_Rbar` (8), `fig3b_Rcrit` (8), `fig3b_bifurcation` (1) | Bifurcation dot `(k_f^c/k_0, R_*/w) = (0.0012020, 2.749)` vs exact Eq. (17) `(0.0012021, 2.750)` |
| **Fig. 3(c)** | 5 | `[800, 2740, 2900, 4160]` | lin / lin; `0.200 px`, `0.350 px` | `fig3c_passive` (10), `fig3c_active` (11) | Finite $V_s/w^3 = 1.8\times 10^6$: active zero drops to `R_bar/w ≈ 22.0` |
| **Fig. 3(d)** | 5 | `[2920, 2740, 4960, 4160]` | log10 / log10; `0.337 px`, `0.428 px` | `fig3d_Rbar` (9), `fig3d_Rcrit` (8) | Finite-system saturation of `R_bar/w` at small $k_f/k_0$ (`33.3` at `1e-6`) |
| **Fig. 4(a)** | 6 | `[780, 1240, 2840, 3380]` | lin / lin; `0.000 px`, `0.000 px` | `fig4a_nullcline_R1` (11), `fig4a_sym_fp` (1) | Open circle `(7.836, 7.837)` matches $V_{\rm tot}\ll V_s$ `7.842` (`0.06%`) vs exact $V_s - V_{\rm tot}$ `7.904`; first x-sample set to `2.2` to clear trajectory arrowhead at `2.0` |
| **Fig. 4(b)** | 6 | `[2960, 1240, 4960, 3380]` | lin / lin; `0.000 px`, `0.000 px` | `fig4b_nullcline_R1` (10), `fig4b_sym_fp` (1) | Filled disk `(6.929, 6.933)` vs Eq. (13)-(14) `6.898` (`0.5%`); nullcline reproduced `10/10` |
| **Fig. 5(a) $V_i$** | 7 | `[800, 1640, 1860, 2580]` | lin / lin; `0.000 px` (2m), `0.350 px` (4m) | `fig5a_V1` (17), `fig5a_V2` (17) | Passive PDE: $V_1/w^3$ grows `2622 -> 3726`, $V_2/w^3$ shrinks `2197 -> 1094` ($\bar{V} \approx 2410\,w^3$) |
| **Fig. 5(a) $\Delta V$** | 7 | `[1860, 1640, 2860, 2580]` | lin / log10; `0.000 px` (2m), `0.218 px` (14m) | `fig5a_dV` (17), `fig5a_dV_fit` (5) | Fit line `(1000, 477.6) -> (8500, 2604.6)` gives $\lambda/k_0 = +2.262\times 10^{-4}$ (Eq. 18 at $\bar{V}=2410$: $+2.897\times 10^{-4}$) |
| **Fig. 5(b)** | 7 | `[2820, 1200, 4940, 2580]` | log10 / log10; `0.360 px` (27m), `0.302 px` (27m) | `fig5b_eq18_line` (6), `fig5b_pde_small` (4), `fig5b_pde_large` (9) | Solid blue line matches Eq. (18) to `0.2%` (`2.33e-4` at `V=3000` vs `2.327e-4`), settling legend typo "Eq. (16)"; PDE markers sit ~24–29% above Eq. (18) |
| **Fig. 5(c) $V_i$** | 7 | `[800, 3180, 1860, 4120]` | lin / lin; `0.000 px` (2m), `0.350 px` (4m) | `fig5c_V1` (17), `fig5c_V2` (17) | Active PDE: $V_1/w^3$ relaxes `5667 -> 4600`, $V_2/w^3$ grows `2591 -> 4300` ($\bar{V} \approx 4450\,w^3$) |
| **Fig. 5(c) $\Delta V$** | 7 | `[1860, 3180, 2860, 4120]` | lin / log10; `0.000 px` (2m), `0.217 px` (14m) | `fig5c_dV` (17), `fig5c_dV_fit` (5) | Fit line `(1000, 2488.4) -> (8500, 351.3)` gives $\lambda/k_0 = -2.610\times 10^{-4}$, revealing 10× caption typo (`5e-3` printed vs `~5e-4` simulated; see `q_digitize_fig5c_caption_rate_typo`) |
| **Fig. 5(d)** | 7 | `[2860, 2640, 4960, 4140]` | log10 / lin; `0.386 px` (24m), `0.167 px` (3m) | `fig5d_eq18_20_line` (11), `fig5d_pde_markers` (11) | Solid blue curve matches Eqs. (18)+(20) (`11/11` reproduced), settling legend typo "Eq. (16)–(17)"; PDE markers cross $\lambda=0$ near $k_b/k_0 \approx 2.6\times 10^{-4}$ |
| **Fig. 6(a)** | 7 | `[5160, 1200, 7300, 2600]` | log10 / log10; `0.342 px` (30m), `0.312 px` (15m) | `fig6a_Nmax` (13), `fig6a_chi1` (7), `fig6a_kfc` (4) | `fig6a_chi1` matches Eq. (21) (`0.6%`), `fig6a_kfc` reads `0.001202` at all 4 heights (`0.00%` vs Eq. 17 `0.0012021`), `fig6a_Nmax` peaks at `27` droplets at `k_f/k_0 = 5.5e-4` |
| **Fig. 6(b)** | 7 | `[7380, 1200, 9300, 2600]` | log10 / log10; `0.326 px` (32m), `0.347 px` (32m) | `fig6b_Bdroplet_upper` (8), `fig6b_Bdroplet_lower` (8), `fig6b_spinodal_upper` (10), `fig6b_spinodal_lower` (9) | Solid B-droplet boundaries match Eq. (17) (`0.1%`) & $k_b=3k_f$ (`0.4%`); dashed spinodal reaches `k_f/k_0 = 0.45 -> k_b/k_0 = 0.520`, confirming **printed** Eq. (F2) (`0.7–1.2%`) and refuting corrected Eq. (F2) (`cutoff 0.125`, 4× lower); upper spinodal runs above $k_b = 3k_f$ for $k_f/k_0 \le 1.5\times 10^{-3}$ |
| **Fig. 7(a)** | 8 | `[5200, 1200, 7240, 2640]` | lin / lin; `0.167 px` (3m), `0.167 px` (3m) | `fig7a_inf` (9), `fig7a_fin` (10), `fig7a_core` (9) | Sharp kink at `R/w = 118.1` on both blue and green curves matches reaction-only Eq. (10) (`0.2–0.3%` med err) vs smooth full Eq. (10)+(12) (`39–45%` med err) |
| **Fig. 7(b)** | 8 | `[7240, 1200, 9320, 2640]` | log10 / log10; `0.246 px` (16m), `0.365 px` (19m) | `fig7b_inf_Rcrit` (4), `fig7b_fin_Rbar` (5), `fig7b_fin_Rcrit` (6), `fig7b_core_Rbar` (11) | `R_crit = 166.9 w` at $\bar{\phi} = 0.011$ confirms $\gamma\beta = w/6$; green curve reaches `R/w = 296.0` ($V = 1.57 V_s$) at $\bar{\phi} = 0.17$, matching `autocatalytic_volume_quad_appB` (`0.5%` med err) and corrected Eq. (B1) (`0.5%` for $\bar{\phi} \ge 0.012$) while printed Eq. (B1) flips $\delta V$ sign at $\bar{\phi} \ge 0.015$ |
| **Fig. 8(a)** | 8 | `[5180, 8240, 7220, 10540]` | lin / lin; `0.167 px` (3m), `0.167 px` (3m) | `fig8a_nullcline_R1` (9), `fig8a_sym_fp` (1) | Open circle `(87.17, 87.18)` matches reaction-only Eq. (10) `87.35` (`0.2%`) vs full Eq. (10) `87.96` (`0.9%`) |
| **Fig. 8(b)** | 8 | `[7300, 8240, 9340, 10540]` | lin / lin; `0.167 px` (3m), `0.167 px` (3m) | `fig8b_nullcline_R1` (8), `fig8b_nullcline_R2` (8), `fig8b_sym_fp` (1) | Filled disk `(93.93, 93.97)` matches reaction-only Eq. (10) `94.12` (`0.2%`) vs full Eq. (10) `90.61` (`3.6%`); full Eq. (10) fails `fig8b_nullcline_R1` (`R_1=115 -> None` vs digitized `49.3`) |
| **Fig. 9(a)** | 9 | `[5180, 1200, 7240, 2620]` | log10 / log10; `0.344 px` (18m), `0.289 px` (20m) | `fig9a_eq24_small` (6), `fig9a_eq24_large` (8), `fig9a_pde_small` (22), `fig9a_pde_large` (20) | Blue Eq. (24) lines (`7.0e-5`, `2.1e-4` at $V=10^4$) match **continuum** $\gamma\beta = w\psi/6 = w/60$ to `0.1–0.2%` (vs **10× error** for $\gamma\beta = w/6$); settles legend typo "Eq. (21)" |
| **Fig. 9(b)** | 9 | `[7220, 1200, 9340, 2620]` | lin / lin; `0.167 px` (3m), `0.150 px` (4m) | `fig9b_Qstab_line` (6), `fig9b_multi_boundary` (8), `fig9b_droplet_boundary` (8) | Red dashed line reads `209.41`, matching Eq. (25) with **capillary** $\gamma\beta = w/6$ (`209.44`, `0.01%`), while continuum `w/60` gives `20.94` (`-90%`) and printed Eq. (24) (`k` instead of `k_s`) gives `2094.4` (`+900%`); black $N_{\max}=1\to 2$ boundary (`187.9–202.3`) sits far above full Eq. (10) (`75–154`, `-56%`) and closely matches `no_kAB` (`184.5–198.0`, `1.6%` med err) |
| **Fig. 10(a)** | 10 | `[800, 1200, 2820, 3180]` | lin / lin; `0.167 px` (3m), `0.000 px` (3m) | `fig10a_nullcline_R1` (7), `fig10a_nullcline_R2` (8), `fig10a_fixed_point` (1) | Filled disk `(116.26, 54.05)` matches reaction-only Eq. (10) `(116.47, 54.16)` (`0.2%`), whereas full Eq. (10)+(12) gives `(109.92, 63.44)` (`+17.4%`, `not reproduced`) |
| **Fig. 10(b)** | 10 | `[2840, 1200, 4980, 2980]` | lin / lin; `0.200 px` (5m), `0.200 px` (4m) | `fig10b_V1` (11), `fig10b_V2` (12), `fig10b_Vtot` (8) | Reaction-only Eq. (10) matches $V_1, V_2, V_{\rm tot}$ to `0.4–0.6%` across `Q_1 = 50..920` with a sharp jump at $Q_1 = Q_{\rm stab} = 209.44$; full Eq. (10)+(12) predicts 2 droplets already at $Q_1 = 150$ inside its `133.5–213.9` bistable window (`20–39%` med err); zero-$\gamma$ $V_1/V_2 = Q_1/Q_2$ deviates by `20–40%` while Eq. (23) matches $V_{\rm tot}$ within `2.8–4.3%` |

---

## Model Comparison Verdicts & Open Questions Settled

### 1. Open Questions Settled Directly by Digitized Figures

1. **`q_figure_legend_equation_typos` — Settled (Confirmed)**
   - `cmp_fig5b_eq18_line` (`close`, 5/6 reproduced, med err `0.2%`): Fig. 5(b)'s solid blue line is **Eq. (18)** (legend misprints "Eq. (16)").
   - `cmp_fig5d_eq18_20_line` (`reproduced`, 11/11): Fig. 5(d)'s solid blue curve is **Eqs. (18)+(20)** (legend misprints "Eq. (16)–(17)").
   - `cmp_fig9a_eq24_small_cont` (`reproduced`, 6/6, med err `0.1%`) & `cmp_fig9a_eq24_large_cont` (`close`, 7/8, med err `0.2%`): Fig. 9(a)'s solid blue lines are **Eq. (24)** (legend misprints "Eq. (21)").

2. **`q_appF_factor_of_2_error` — Settled (Confirmed that Fig. 6b plotted printed Eq. F2)**
   - `cmp_fig6b_spinodal_upper_printed` (`close`, 8/10 rep, med err `0.7%`) and `cmp_fig6b_spinodal_lower_printed` (`close`, 7/9 rep, med err `1.2%`) closely reproduce Fig. 6(b)'s dashed spinodal boundary up to $k_f/k_0 = 0.45$ ($k_b/k_0 = 0.520$, apex at $0.5, 0.5$).
   - `cmp_fig6b_spinodal_upper_corrected` (`not reproduced`, med err `17.5%`, no solution for $k_f/k_0 > 0.125$) and `cmp_fig6b_spinodal_lower_corrected` (`not reproduced`, med err `43.8%`) prove quantitatively that Fig. 6(b) was generated with the erroneous printed Eq. (F2) (cutoff $k_f + k_b < k_0$) rather than the exact linearization of Eq. (3b) (cutoff $k_f + k_b < k_0/4$, verified by `binary_cahn_hilliard_1d_pde`).

3. **`q_fig6b_spinodal_vs_3kf_overlap` — Settled (Confirmed)**
   - At $k_f/k_0 \in \{6\times 10^{-4}, 10^{-3}, 1.5\times 10^{-3}\}$, `fig6b_spinodal_upper_data` reads $k_b/k_0 = \{0.00215, 0.00355, 0.00527\}$, sitting strictly above `fig6b_Bdroplet_lower_data` ($k_b = 3k_f = \{0.00180, 0.00301, 0.00450\}$) inside the blue B-rich droplet wedge, matching the analytical low-rate asymptote $k_b/k_f = 2 + \sqrt{3} \approx 3.732$.

4. **`q_sec4_length_scale_w_inconsistency` — Settled (Confirmed 10× split within Sec. IV)**
   - **Fig. 9(a)** (`cmp_fig9a_eq24_small_cont`: `reproduced`, `0.1%`; `cmp_fig9a_eq24_large_cont`: `close`, `0.2%`) uses continuum $\gamma\beta = w\psi/6 = w/60$, whereas capillary $\gamma\beta = w/6$ fails by exactly **10×** (`cmp_fig9a_eq24_small_cap`, `cmp_fig9a_eq24_large_cap`: both `not reproduced`, `+900%`).
   - **Fig. 9(b)** (`cmp_fig9b_Qstab_cap`: `reproduced`, 6/6, `0.01%`) and **Figs. 7, 8, 10** (`cmp_fig7a_inf`, `cmp_fig7b_inf_Rcrit`, `cmp_fig8a/b`, `cmp_fig10a/b`) use capillary $\gamma\beta = w/6$, whereas continuum $\gamma\beta = w/60$ fails by **10×** (`cmp_fig9b_Qstab_cont`: `not reproduced`, `-90%`).

5. **`q_eq24_missing_subscript_s` — Settled (Confirmed on Fig. 9b)**
   - Solving Eq. (24) $\lambda = 0$ with the **corrected** $k_s = k\psi$ yields Eq. (25) $Q_{\rm stab}/(w^3 k_0) = 209.44$ (`cmp_fig9b_Qstab_cap`: `reproduced`, `0.01%` vs digitized `209.41`), whereas solving **printed** Eq. (24) with $k$ instead of $k_s$ yields $2094.4$ (`cmp_fig9b_Qstab_eq24_printed`: `not reproduced`, `+900%`).

6. **`q_eqB1_surface_tension_term_error` — Settled (Confirmed on Fig. 7b)**
   - On `fig7b_core_Rbar_data` ($\bar{\phi} \ge 0.012$), `autocatalytic_volume_quad_appB` (`cmp_fig7b_core_Rbar_quad`: `reproduced`, `0.5%` med err) and **corrected** Eq. (B1) (`cmp_fig7b_core_Rbar_eqB1_corrected`: `0.5%` med err for $\bar{\phi} \ge 0.012$, positive $\delta V$ above Eq. 23) track the plotted green curve; **printed** Eq. (B1) (`cmp_fig7b_core_Rbar_eqB1_printed`) erroneously flips the sign of $\delta V$ below Eq. (23) for $\bar{\phi} \ge 0.015$ (e.g., `93.61` vs Eq. 23 `93.74` and digitized `94.26` at $\bar{\phi} = 0.015$). Both Eq. (B1) variants are undefined for $\bar{\phi} \le k_b/k = 0.010$ where $V_0 \le 0$.

7. **`q_audit_sec4_autocat_bg_bistability` & New `q_digitize_sec4_reaction_only_bg` — Settled**
   - Full Eq. (10)+(12) (`autocatalytic_emulsion_solver`) disagrees systematically with the plotted Sec. IV figures: `cmp_fig8b_nullcline_R1_full` (`not reproduced`), `cmp_fig10a_nullcline_R1_full` (`not reproduced`), `cmp_fig10a_fp_full` (`not reproduced`, `63.44` vs `54.05`, `+17.4%`), `cmp_fig7b_core_Rbar_full` (`not reproduced`, undefined for $\bar{\phi} \ge 0.11$ where $V > V_s$), `cmp_fig9b_multi_full` (`not reproduced`, `75–154` vs `188–202`, `-55.9%`), and baseline median errors of `39–45%` on Fig. 7(a) and `20–39%` on Fig. 10(b) (whose BoxGraph verdicts widen to `reproduced` only via `_carry_rounding` on 1-sig-fig `0.02` and `2e-4` inputs).
   - By contrast, the reaction-only, dilute-volume closure ($\phi_0^B = \max(0, \bar{\phi} - \psi V_{\rm tot}/V_s - k_b/k)$) reproduces **all** of Figs. 7(a) (including the sharp kink at $R/w = 118.1$), 7(b) (including $V > V_s$ at $\bar{\phi} = 0.17$), 8(a), 8(b), 10(a) (`(116.47, 54.16)` vs `(116.26, 54.05)`), and 10(b) (`0.2–0.6%` median error throughout) and eliminates full Eq. (10)'s `133.5–213.9` bistability window.
   - Similarly, in passive **Fig. 4(a)**, the dilute-volume closure $\phi_0^B = \bar{\phi}^B - V_{\rm tot}/V_s$ reproduces `fig4a_sym_fp` (`7.842` vs `7.837`, `0.06%`) and `fig4a_nullcline_R1` (`11/11` reproduced, `0.1%` med err), whereas exact $(V_s - V_{\rm tot})$ shifts $R_{\rm sym}$ to `7.904` (`+0.9%`).

8. **New `q_digitize_fig5c_caption_rate_typo` — Discovered via Digitization**
   - Fig. 5(c)'s caption prints $k_b/k_0 = 5\times 10^{-3}, k_f/k_0 = 5\times 10^{-4}$, but `fig5c_dV_fit_data` yields $\lambda/k_0 = -2.61\times 10^{-4}$ (10× slower than Eq. 18's $-3.16\times 10^{-3}$ and `run_pde_axi_fig5c_active`'s $-2.50\times 10^{-3}$) and `fig5c_V1/V2` converge toward $V_1 \approx V_2 \approx 4500\,w^3$ (exceeding Eq. 16's maximum single-droplet volume $1687\,w^3$ at $k_b/k_0 = 5\times 10^{-3}$). Both match 10× smaller rates ($k_b/k_0 \approx 5\times 10^{-4}, k_f/k_0 \approx 5\times 10^{-5}$), consistent with Fig. 5(d)'s x-axis range.

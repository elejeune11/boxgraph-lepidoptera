# Figure Inventory for `claessens2006`

All figures in `corpus/actin_bundling/2006_Claessens.pdf` (PDF pp. 1–19: manuscript, Methods, References, Captions, Figs. 1–3; PDF p. 20: Supplementary Information NM06040751A, SI Fig. 1a–b) were inspected and logged below. Every quantitative panel (Fig. 3A–D and SI Fig. 1a–b) was digitized with self-contained calibrated Python reading scripts (`figures/*.py`, `mode: code`), checked via `boxgraph figure check`, visually verified against contact sheets of the rendered overlays, and saved in `readings.json`.

**Summary**: 23 calibrated figure readings across 6 quantitative panels — 105 curve points across 11 solid/dashed model curves + 188 measured data points across 12 experimental series = **293 total digitized points**. No plotted data series carry error bars.

## Figures Digitized (in `figures.yaml` and `figures/`)

- **Fig. 3A (p. 19, caption p. 16) — Plastin-crosslinked bundles ($\kappa_B$ vs $n$, log–log, $n \in [1, 40]$, $\kappa_B \in [0.04, 80]\ \text{pN}\cdot\mu\text{m}^2$)**:
  - Rendered at 900 dpi (crop `[900, 1830, 3750, 4140]`; native embedded 150 ppi JPEG, 6 rendered px per native px). Calibrated from all 13 x-ticks ($n = 1, 2\text{--}10, 20, 30, 40$; RMS misfit $1.85\text{ px} = 0.31$ native px) and 59 left/right y-ticks ($0.05\text{--}70\ \text{pN}\cdot\mu\text{m}^2$; RMS misfit $2.10\text{ px} = 0.35$ native px).
  - `fig3a_model_plastin_05` (`fig3a_plastin_05_curve`, 9 pts at $n = 2, 2.5, 10, 15, 20, 25, 30, 35, 38$; $\kappa_B = 0.153\text{--}4.24\ \text{pN}\cdot\mu\text{m}^2$): red solid model curve for $r_{\text{ABP:a}} = 1{:}2$ ($0.5$).
  - `fig3a_model_plastin_002` (`fig3a_plastin_002_curve`, 10 pts at $n = 2, 6, 8, 10, 15, 20, 25, 30, 35, 38$; $\kappa_B = 0.142\text{--}2.73\ \text{pN}\cdot\mu\text{m}^2$): dark-blue solid model curve for $r_{\text{ABP:a}} = 1{:}50$ ($0.02$), lying closely atop the lower decoupled dashed line (`reproduced` by variant $K = S_2$, `close` under printed $K = S_2/N$).
  - `fig3a_asymptote_coupled` (`fig3a_coupled_asymptote_curve`, 9 pts at $n = 2, 3, 5, 8, 10, 15, 20, 25, 30$; $\kappa_B = 0.284\text{--}62.9\ \text{pN}\cdot\mu\text{m}^2$): upper black dashed line $\kappa_B = n^2 \kappa_a$ ($t = 0$). Reproduces $\kappa_a = 7 \times 10^{-26}\ \text{N}\cdot\text{m}^2$ (p. 5; `reproduced`), whereas $7.3 \times 10^{-26}\ \text{N}\cdot\text{m}^2$ (p. 12) is ~4% high (`close`).
  - `fig3a_data_plastin_05` (`fig3a_plastin_05_meas`, 7 visible red filled circles at $n = 1.19\text{--}11.8$, $\kappa_B = 0.0653\text{--}1.21\ \text{pN}\cdot\mu\text{m}^2$): measured $r_{\text{ABP:a}} = 1{:}2$ (caption prints swapped circle colors; a few additional red dots are heavily occluded behind white-filled blue open circles in the $n \approx 6\text{--}9$ cluster and cannot be uniquely separated).
  - `fig3a_data_plastin_002` (`fig3a_plastin_002_meas`, 14 blue open circles at $n = 2.94\text{--}12.3$, $\kappa_B = 0.207\text{--}1.09\ \text{pN}\cdot\mu\text{m}^2$): measured $r_{\text{ABP:a}} = 1{:}50$.

- **Fig. 3B (p. 19, caption p. 16) — Fascin-crosslinked bundles ($\kappa_B$ vs $n$, log–log)**:
  - Rendered at 900 dpi (crop `[3810, 1830, 6690, 4140]`). Calibrated from 13 x-ticks (RMS $1.83\text{ px}$) and 59 y-ticks (RMS $2.40\text{ px}$).
  - `fig3b_model_fascin_05` (`fig3b_fascin_05_curve`, 9 pts at $n = 2, 3, 5, 8, 10, 15, 25, 35, 38$; $\kappa_B = 0.276\text{--}33.3\ \text{pN}\cdot\mu\text{m}^2$): top red solid curve ($r = 0.5$).
  - `fig3b_model_fascin_02` (`fig3b_fascin_02_curve`, 10 pts at $n = 2, 3, 5, 9, 10, 15, 20, 25, 30, 34$; $\kappa_B = 0.230\text{--}15.2\ \text{pN}\cdot\mu\text{m}^2$): second dark-blue solid curve. Comparison settles `q_caption_vs_si_concentration_ratios`: matches SI ratio $r = 0.2$ (within $-11.4\%$ to $+2.3\%$ for $n \ge 9$ under $K = S_2$), whereas caption ratio $1{:}10$ ($0.1$) underpredicts by $-36\%$ to $-41\%$.
  - `fig3b_model_fascin_005` (`fig3b_fascin_005_curve`, 10 pts at $n = 2, 3, 5, 8, 10, 15, 20, 30, 35, 38$; $\kappa_B = 0.189\text{--}8.51\ \text{pN}\cdot\mu\text{m}^2$): third green solid curve ($r = 0.05$).
  - `fig3b_model_fascin_002` (`fig3b_fascin_002_curve`, 10 pts at $n = 3, 5, 8, 12, 15, 20, 25, 30, 35, 38$; $\kappa_B = 0.251\text{--}4.48\ \text{pN}\cdot\mu\text{m}^2$): lowest orange solid curve ($r = 0.02$; `close` under $K = S_2$, within $5.8\%$ for $n \ge 9$, vs `+64%` overshoot under printed $K = S_2/N$).
  - `fig3b_data_fascin_05` (`fig3b_fascin_05_meas`, 15 red filled circles at $n = 1.99\text{--}30.5$, $\kappa_B = 0.148\text{--}27.3\ \text{pN}\cdot\mu\text{m}^2$).
  - `fig3b_data_fascin_02` (`fig3b_fascin_02_meas`, 14 blue open circles at $n = 3.25\text{--}36.1$, $\kappa_B = 0.370\text{--}17.6\ \text{pN}\cdot\mu\text{m}^2$).
  - `fig3b_data_fascin_005` (`fig3b_fascin_005_meas`, 10 green open squares at $n = 5.02\text{--}32.8$, $\kappa_B = 0.855\text{--}7.27\ \text{pN}\cdot\mu\text{m}^2$; confirms `q_fig3b_caption_square_symbols_swapped`: green `□` is $r = 0.05$, above orange `■`).
  - `fig3b_data_fascin_002` (`fig3b_fascin_002_meas`, 13 orange filled squares at $n = 2.66\text{--}32.2$, $\kappa_B = 0.142\text{--}3.11\ \text{pN}\cdot\mu\text{m}^2$; confirms orange `■` is $r = 0.02$).
  - Dashed asymptotes omitted on Fig. 3B as identical reference lines are read cleanly on Fig. 3A (upper) and Fig. 3D (lower).

- **Fig. 3C (p. 19, caption p. 16) — $\alpha$-actinin-crosslinked bundles ($\kappa_B$ vs $n$, log–log)**:
  - Rendered at 900 dpi (crop `[900, 4605, 3750, 6945]`). Calibrated from 13 x-ticks (RMS $2.55\text{ px}$) and 57 y-ticks (RMS $2.22\text{ px}$).
  - `fig3c_model_aact_01` (`fig3c_aact_01_curve`, 9 pts at $n = 3, 5, 8, 15, 20, 25, 30, 35, 38$; $\kappa_B = 0.713\text{--}31.8\ \text{pN}\cdot\mu\text{m}^2$): red solid curve ($r = 0.1$).
  - `fig3c_model_aact_002` (`fig3c_aact_002_curve`, 10 pts at $n = 3, 4, 6, 12, 15, 20, 25, 30, 35, 38$; $\kappa_B = 0.288\text{--}7.16\ \text{pN}\cdot\mu\text{m}^2$): dark-blue solid curve. Settles `q_caption_vs_si_concentration_ratios`: consistent with caption ratio $1{:}50$ ($r = 0.02$), whereas SI table ratio $r = 0.05$ overshoots by $+81\%$ to $+205\%$ under $K = S_2$ and $+120\%$ to $+586\%$ under printed $K = S_2/N$.
  - `fig3c_data_aact_01` (`fig3c_aact_01_meas`, 41 red filled circles at $n = 1.00\text{--}28.6$, $\kappa_B = 0.0604\text{--}14.4\ \text{pN}\cdot\mu\text{m}^2$).
  - `fig3c_data_aact_002` (`fig3c_aact_002_meas`, 15 blue open circles at $n = 1.05\text{--}10.1$, $\kappa_B = 0.0738\text{--}1.09\ \text{pN}\cdot\mu\text{m}^2$).
  - Dashed asymptotes omitted on Fig. 3C as identical reference lines are read cleanly on Fig. 3A (upper) and Fig. 3D (lower).

- **Fig. 3D (p. 19, caption p. 16) — PEG 6 kDa depletion bundles ($\kappa_B$ vs $n$, log–log)**:
  - Rendered at 900 dpi (crop `[3810, 4605, 6690, 6945]`). Calibrated from 13 x-ticks (RMS $2.25\text{ px}$) and 59 y-ticks (RMS $2.41\text{ px}$).
  - `fig3d_model_peg` (`fig3d_peg_curve`, 9 pts at $n = 2, 3.8, 10, 13, 18.5, 27, 30, 32, 38$; $\kappa_B = 0.254\text{--}54.6\ \text{pN}\cdot\mu\text{m}^2$): red solid PEG model curve ($t = 0$, $k_\parallel / \delta = 10^3\ \text{N/m}^2$, $L_B = 21.2 + 1.3 n\ \mu\text{m}$, independent of unstated $c_a$). **Decisive test of $K$**: printed $K = S_2/N$ (`cmp_fig3d_peg_curve_printed`) is `not reproduced` ($+10.7\%$ to $+44.5\%$ overshoot for $n = 13\text{--}38$), whereas variant $K = S_2$ (`cmp_fig3d_peg_curve_KnoN`, overall verdict `close`) lands inside the 95% reading interval of all 7 points from $n = 10$ to $n = 38$ ($+1.4\%$ to $+3.8\%$), differing only at $n = 2, 3.8$ ($+13\%$ to $+14\%$) where $N < 1$.
  - `fig3d_asymptote_decoupled` (`fig3d_decoupled_asymptote_curve`, 10 pts at $n = 2, 3, 5, 8, 10, 15, 20, 25, 30, 35$; $\kappa_B = 0.137\text{--}2.46\ \text{pN}\cdot\mu\text{m}^2$): lower black dashed line $\kappa_B = n \kappa_a$, reproducing $\kappa_a = 7 \times 10^{-26}\ \text{N}\cdot\text{m}^2$ (p. 5; `reproduced`) vs $7.3 \times 10^{-26}\ \text{N}\cdot\text{m}^2$ (`close`).
  - `fig3d_data_peg_4pct` (`fig3d_peg_4pct_meas`, 26 red filled circles at $n = 1.59\text{--}33.5$, $\kappa_B = 0.0755\text{--}44.6\ \text{pN}\cdot\mu\text{m}^2$): 4% w/w PEG 6 kDa.
  - `fig3d_data_peg_2pct` (`fig3d_peg_2pct_meas`, 15 blue open circles at $n = 2.71\text{--}28.0$, $\kappa_B = 0.443\text{--}38.1\ \text{pN}\cdot\mu\text{m}^2$): 2% w/w PEG 6 kDa.

- **SI Fig. 1a (p. 20, printed SI 1) — Dynamic equilibration of fuzzy diameter $D$ vs time averaged (fascin:actin 1:2, $R_{\text{ring}} = 8.0\ \mu\text{m}$)**:
  - Rendered at 1200 dpi (crop `[1220, 5200, 4440, 7680]`). Calibrated linearly from all 9 x-ticks ($t = 0\text{--}8\ \text{s}$; RMS $1.71\text{ px}$) and 7 y-ticks ($D = -0.02\text{ to }0.10\ \mu\text{m}$; RMS $1.63\text{ px}$).
  - `si_fig1a_D_vs_t` (`si_fig1a_D_vs_t_meas`, 10 red points at $t = 0.112\text{--}7.03\ \text{s}$, exact multiples $1, 3, 5, 10, 15, 20, 30, 40, 50, 60$ of the 117 ms frame interval; $D = 0.0000\text{--}0.0814\ \mu\text{m}$, saturating cleanly at $D = 0.0814 \pm 0.0003\ \mu\text{m}$ for $t \ge 3.53\ \text{s}$). Implies $\kappa_B = 50.9\ \text{pN}\cdot\mu\text{m}^2$, confirming `q_si_fig1a_stiffness_vs_fig3b` against Fig. 3B fascin 1:2 ($\max \kappa_B = 27.3\ \text{pN}\cdot\mu\text{m}^2$). No analytical $D(t)$ model is given in the paper, so recorded as an uncompared measured curve for downstream validation.

- **SI Fig. 1b (p. 20, printed SI 1) — Dynamic equilibration of fuzzy diameter $D$ vs time averaged (fascin:actin 1:50, $R_{\text{ring}} = 5.2\ \mu\text{m}$)**:
  - Rendered at 1200 dpi (crop `[5100, 5200, 8320, 7680]`). Calibrated linearly from all 9 x-ticks ($t = 0\text{--}8\ \text{s}$; RMS $1.74\text{ px}$) and 8 y-ticks ($D = -0.05\text{ to }0.30\ \mu\text{m}$; RMS $1.98\text{ px}$).
  - `si_fig1b_D_vs_t` (`si_fig1b_D_vs_t_meas`, 8 red points at $t = 0.131\text{--}7.02\ \text{s}$, multiples $1, 5, 10, 20, 30, 40, 50, 60$ of 117 ms; $D = -0.0009\text{--}0.2571\ \mu\text{m}$, saturating cleanly at $D = 0.2571 \pm 0.0008\ \mu\text{m}$ for $t \ge 4.68\ \text{s}$). Implies $\kappa_B = 1.40\ \text{pN}\cdot\mu\text{m}^2$, consistent with Fig. 3B fascin 1:50 ($0.142\text{--}3.11\ \text{pN}\cdot\mu\text{m}^2$).

## Figures Left Out (and Why)

- **Fig. 1 (p. 17, caption p. 16)**: Schematic diagram of the discrete ABP-crosslinked F-actin bundle model showing filament diameter $d$, inter-filament spacing $t$, axial spacing $\delta$, shear spring $k_\parallel$, and decoupled ($\alpha \ll 1$) vs fully coupled ($\alpha \gg n$) bending deformation. Qualitative illustration with no numeric axes or plotted data.
- **Fig. 2a, 2b (p. 18, caption p. 16)**: Fluorescence microscopy images of F-actin inside water-in-dodecane emulsion droplets without ABPs (2a) and organized into a closed bundle ring of radius $R_r$ and fuzzy diameter $D$ inside a droplet of radius $R_d$ (2b), with $10\ \mu\text{m}$ scale bars. Micrographs with no numeric plot axes or data series.

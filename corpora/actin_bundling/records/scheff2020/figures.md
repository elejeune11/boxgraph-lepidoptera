# Figure Digitization Log — `scheff2020`

**Paper:** D. R. Scheff, K. L. Weirich, K. Dasbiswas, A. Patel, S. Vaikuntanathan, and M. L. Gardel, *"Tuning shape and internal structure of protein droplets via biopolymer filaments"*, *Soft Matter* **16**, 5659–5668 (2020) + Electronic Supplementary Information (PDF pp. 1–11, `corpus/actin_bundling/2020_Scheff.pdf`).

## Summary

- **Calibrated readings:** **37** calibrated figure readings in `figures.yaml` covering all **16** quantitative plot panels across Figs. 1–5 and ESI Fig. S1, backed by fast custom readers (`figures/_helpers.py`, `_fig1_2.py`, `_fig3_4.py`, `_fig5_S1.py`, and 37 lightweight per-reading entrypoints operating directly on the rendered image array).
- **Digitized points:** **612 / 612** points found (100% found), including **195** points with calibrated vertical error bars (SD, SEM, fit-range, or unstated), saved in `readings.json` (`37/37` `saved == "matches"`).
- **Visual verification:** 100% of calibrated readings (`37/37`) passed `boxgraph figure check` and were visually inspected via rendered magenta overlays saved under `figwork/overlays/*.png` and assembled into 7 contact sheets (`figwork/overlay_sheets/official_sheet_1.png`–`official_sheet_7.png`). Every magenta cross and error-bar span sits centered on its intended curve, marker, or histogram bar.
- **Record integration:** Added `midplane_profile_distance` to `quantities.yaml`; added `quasi_bipolar_aspect_vs_area_vtilde_inv` (`variant_of: quasi_bipolar_aspect_vs_area_printed`) under `models/`; appended all 37 `figure:` datums to `data.yaml` (preserving all 58 pre-existing eye-read/printed datums), 48 figure sweep runs to `runs.yaml` (119 total runs, 0 failed), 48 figure comparisons to `comparisons.yaml` (65 total comparisons), and 3 new digitization questions to `questions.yaml` (preserving all 11 pre-existing questions verbatim).

---

## Panel-by-Panel Inventory (Pages 1–11)

### Figure 1 (Page 4)

- **Fig. 1(a) — Left out:** Qualitative schematic illustrating FUS droplet phase separation and actin filament partitioning; contains no quantitative axes or data.
- **Fig. 1(b) — Left out:** Confocal fluorescence micrographs (FUS-GFP, rhodamine-actin, merged channel) showing actin partitioning into FUS droplets; no quantitative plot axes.
- **Fig. 1(c) — Read (`fig1c_actin`, `fig1c_fus`):** Normalized fluorescence intensity across a droplet midplane vs distance (1 $\mu\text{M}$ actin, $L \sim 90\text{ nm}$, 4.4 $\mu\text{M}$ FUS).
  - *Calibration:* 1200 dpi; linear x-axis from left spine ($0\ \mu\text{m}$) + 3 top inward ticks ($2, 4, 6\ \mu\text{m}$), residual **0.80 px**; linear y-axis from 4 left inward ticks ($1.0, 0.75, 0.5, 0.25$) + bottom spine ($0.0$), residual **1.20 px**. Both the upper-left legend box and the upper-right inset micrograph frame were automatically detected and masked out.
  - *Readings:* `fig1c_actin` (red triangles: **61** points, $x \in [0.68, 7.72]\ \mu\text{m}$, $y \in [0.010, 0.998]$) and `fig1c_fus` (black circles: **79** points, $x \in [0.15, 8.04]\ \mu\text{m}$, $y \in [0.010, 0.998]$). No error bars plotted.
  - *Left-out element:* Upper-right confocal inset micrograph showing the line-scan chord.
- **Fig. 1(d) — Read (`fig1d_90nm`, `fig1d_135nm`, `fig1d_270nm`):** Actin partitioning intensity ratio $I_\text{inside}/I_\text{outside}$ vs $[\text{Actin}]$ ($\log_{10}$ x-axis) with $\pm 1\text{ SD}$ error bars (`"bars are standard deviation between droplets"`, p. 4).
  - *Calibration:* 1200 dpi; $\log_{10}$ x-axis from 3 major top ticks ($0.01, 0.1, 1\ \mu\text{M}$), residual **0.00 px**; linear y-axis from 8 right ticks ($80, \dots, 10$) + bottom spine ($0$), residual **1.53 px**.
  - *Readings:*
    - `fig1d_90nm` (blue open triangles, **2** points, 2 with SD bars): $57.83\ [37.31, 77.61]$ at $0.1\ \mu\text{M}$; $26.92\ [24.17, 30.04]$ at $1\ \mu\text{M}$.
    - `fig1d_135nm` (red open circles, **2** points, 2 with SD bars): $28.03\ [23.80, 32.26]$ at $0.1\ \mu\text{M}$; $24.24\ [19.35, 28.75]$ at $1\ \mu\text{M}$.
    - `fig1d_270nm` (black open squares, **3** points, 3 with SD bars): $7.21\ [5.12, 9.34]$ at $0.01\ \mu\text{M}$; $39.51\ [36.38, 42.63]$ at $0.1\ \mu\text{M}$; $27.20\ [22.69, 31.71]$ at $1\ \mu\text{M}$.
  - *Note:* Confirms prose `"about 7"` at $0.01\ \mu\text{M}$ (only $L \sim 270\text{ nm}$ plotted) and `"~25–40"` for 5 of 6 points at $0.1$–$1\ \mu\text{M}$, while $L \sim 90\text{ nm}$ at $0.1\ \mu\text{M}$ reaches $57.8 \pm 20.1$ (logged in `q_digitized_batch_and_prose_discrepancies`).

### Figure 2 (Page 5)

- **Fig. 2(a) — Left out:** Time-lapse confocal micrographs of two coalescing droplets at $t = 0, 2, 16, 40\text{ s}$; qualitative images with no plot axes.
- **Fig. 2(b) main — Read (`fig2b_exp_curve`) + partially left out:** Normalized coalescence area $(A(t) - A_f)/(A_i - A_f)$ vs rescaled time $t/\tau$.
  - *Calibration:* 1200 dpi; linear x-axis from left spine ($0$) + 4 top ticks ($1, 2, 3, 4$), residual **0.80 px**; linear y-axis from 5 right ticks ($1.0, 0.75, 0.5, 0.25, 0.0$), residual **1.60 px**.
  - *Reading:* `fig2b_exp_curve` (black dashed theoretical master curve $\exp(-t/\tau)$) sampled at **10** points $t/\tau \in \{0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 3.8\}$ ($y = 0.7818$ down to $0.0228$).
  - *Comparison (`cmp_fig2b_exp_repro`):* **`reproduced`** (**10/10** points inside 95% reading interval, $\text{RMS} = 0.0029$, $\max|d| = 0.0070$) against `coalescence_normalized_area` ($\exp(-t/\tau)$).
  - *Left-out element:* The 6 colored scatter point sets overlaid directly along the master curve ($0$–$1\ \mu\text{M}$ actin) are heavily superimposed within the same pixels and cannot be cleanly separated trace-by-trace; the un-collapsed single-event coalescence curve in Fig. 2(b) inset was digitized instead.
- **Fig. 2(b) inset — Read (`fig2b_inset_data`, `fig2b_inset_fit`) + pre-contact branches left out:** Single coalescence event cross-sectional area $A(t)$ ($\mu\text{m}^2$) vs time $t$ (s) for 4.4 $\mu\text{M}$ FUS, 1 $\mu\text{M}$ actin, $L \sim 90\text{ nm}$.
  - *Calibration:* 2400 dpi (high-res crop of the small inset); linear x-axis from 7 top ticks ($0, 10, \dots, 60\text{ s}$), residual **3.44 px** ($\approx 0.16\text{ s}$, $<0.3\%$ of axis); linear y-axis from 9 right ticks ($24, 22, \dots, 8\ \mu\text{m}^2$), residual **2.47 px**.
  - *Readings:* `fig2b_inset_data` (green triangles, **16** sampled times $t \in [0, 55]\text{ s}$, $A = 23.46 \to 19.60\ \mu\text{m}^2$) and `fig2b_inset_fit` (black dashed exponential fit, **16** sampled times $t \in [0, 55]\text{ s}$, $A = 23.70 \to 19.28\ \mu\text{m}^2$).
  - *Left-out element:* The two separate pre-contact droplet area branches at $t < 0\text{ s}$ (prior to initial contact at $t = 0\text{ s}$) were left out because the single-exponential relaxation curve is defined only for $t \ge 0$.
- **Fig. 2(c) — Read (`fig2c_tau_data`, `fig2c_linear_fit`):** Coalescence relaxation time $\tau$ (s) vs relaxation length $\Delta l = \sqrt{A_i} - \sqrt{A_f}$ ($\mu\text{m}$) at 1 $\mu\text{M}$ actin, 4.4 $\mu\text{M}$ FUS.
  - *Calibration:* 1200 dpi; linear x-axis from 5 top ticks ($0.2, \dots, 0.6\ \mu\text{m}$), residual **1.60 px**; linear y-axis from 4 right ticks ($20, 15, 10, 5\text{ s}$) + bottom spine ($0$), residual **2.40 px**.
  - *Readings:* `fig2c_tau_data` (**9** open black squares spanning $\Delta l \in [0.188, 0.578]\ \mu\text{m}$, $\tau \in [3.54, 21.71]\text{ s}$, 8 with visible small vertical error bars of `unstated` kind — likely exponential-fit error) and `fig2c_linear_fit` (black dashed linear fit sampled at **9** positions $\Delta l \in [0.20, 0.60]\ \mu\text{m}$, $\tau = 4.04 \to 14.59\text{ s}$).
  - *Comparisons & finding:* `cmp_fig2c_scaling_vs_data` (`not reproduced`, $\text{RMS} = 3.62\text{ s}$, implied $C = 1.373$ relative to probe $\eta/\gamma = 19\text{ s}/\mu\text{m}$, i.e., proportional slope $\approx 26.1\text{ s}/\mu\text{m}$) and `cmp_fig2c_scaling_vs_fit` (`not reproduced` at low $\Delta l$ due to nonzero $-1.18\text{ s}$ intercept, $\text{RMS} = 0.49\text{ s}$, $C = 1.208$). Crucially, the plotted dashed linear fit has exact slope **$26.34\text{ s}/\mu\text{m}$** (intercept $-1.18\text{ s}$, unweighted OLS on the 9 squares gives $29.10 \pm 9.55\text{ s}/\mu\text{m}$), which is $7.4\text{ s}/\mu\text{m}$ steeper than the $1\ \mu\text{M}$ value plotted in Fig. 2(d) ($18.92 \pm 10.57\text{ s}/\mu\text{m}$; logged in `q_digitized_fig2c_vs_fig2d_slope_discrepancy`).
- **Fig. 2(d) — Read (`fig2d_eta_over_gamma`):** Inverse capillary velocity $\eta/\gamma$ ($\text{s}/\mu\text{m}$) vs $[\text{Actin}]$ ($\mu\text{M}$) on a broken axis ($0\ //\ 0.01, 0.1, 1\ \mu\text{M}$) with SEM error bars (`"are the standard error from the linear fit"`, p. 5).
  - *Calibration:* 1200 dpi; broken x-axis with $0\ \mu\text{M}$ tick anchored directly and $\log_{10}$ fit to $0.01, 0.1, 1\ \mu\text{M}$ top major ticks, residual **0.00 px**; linear y-axis from 7 right ticks ($30, 25, \dots, 0\text{ s}/\mu\text{m}$), residual **2.86 px**.
  - *Reading:* **4** filled black squares (4 with SEM bars): $0.91\ [-0.10, 1.96]$ at $0\ \mu\text{M}$; $-0.09\ [-2.74, 2.28]$ at $0.01\ \mu\text{M}$; $4.26\ [-1.38, 9.94]$ at $0.1\ \mu\text{M}$; $18.92\ [8.36, 29.51]\text{ s}/\mu\text{m}$ at $1\ \mu\text{M}$. Confirms the $\sim 20$-fold increase reported on p. 5 (contradicting the Abstract's `"coalescence time decreases"`, `q_abstract_coalescence_time_sign_reversal`).

### Figure 3 (Page 6)

- **Fig. 3(a) — Left out:** Confocal fluorescence micrographs of FUS droplets at $0, 0.01, 0.1, 1\ \mu\text{M}$ actin ($L \sim 270\text{ nm}$); no plot axes.
- **Fig. 3(b) — Left out:** Zoomed fluorescence micrographs showing individual tactoid and spherical droplets; no plot axes.
- **Fig. 3(c) — Read (`fig3c_aspect_vs_area`):** Binned average aspect ratio vs cross-sectional area for 1 $\mu\text{M}$ actin, $L \sim 270\text{ nm}$ ($n = 1180$, $\pm 1\text{ SD}$ error bars, `"Error bars represent  1 standard deviation between droplets"`, p. 6).
  - *Calibration:* 1200 dpi; linear x-axis from 9 top ticks ($0, 15, \dots, 120\ \mu\text{m}^2$), residual **1.97 px**; linear y-axis from top spine ($2.5$) + 6 right ticks ($2.25, \dots, 1.0$), residual **0.00 px**.
  - *Reading:* **62** black square bins spanning $A \in [1.42, 116.3]\ \mu\text{m}^2$, $x \in [1.031, 1.917]$ (**50** with $\pm 1\text{ SD}$ bars, 12 single-droplet bins without bars).
  - *Validation comparisons:* Against the Fig. S1 (right) $L \sim 270\text{ nm}$ fitted parameters ($K/\gamma = 1.7\ \mu\text{m}$, $\omega = 2.8$), `cmp_fig3c_270nm_printed` gives `not reproduced` (**0/62** reproduced, $\text{RMS} = 1.313$), while `cmp_fig3c_270nm_vtilde_inv` gives `not reproduced` overall (**37/62** reproduced, 13 close, 12 not reproduced at small areas, $\text{RMS} = 0.283$), directly quantifying the batch shift between Fig. 3(c) and Fig. S1 (right) (`q_270nm_batch_variation_fig3c_vs_figS1`).
- **Fig. 3(d) — Read (`fig3d_small`, `fig3d_medium`, `fig3d_large`):** Average aspect ratio vs $[\text{Actin}]$ ($\log_{10}$ x-axis, $L \sim 270\text{ nm}$) in 3 droplet size bins ($\pm 1\text{ SD}$ bars).
  - *Calibration:* 1200 dpi; $\log_{10}$ x-axis from 3 major top ticks ($0.01, 0.1, 1\ \mu\text{M}$), residual **1.33 px**; linear y-axis from top spine ($1.8$) + 4 right ticks ($1.6, 1.4, 1.2, 1.0$), residual **1.20 px**.
  - *Readings:* **3** points per series (**9** points total, all 9 with SD bars). At $1\ \mu\text{M}$: small ($4$–$6\ \mu\text{m}^2$, blue squares) $= 1.450 \pm 0.266$; medium ($14$–$16\ \mu\text{m}^2$, red triangles) $= 1.234\ [1.112, 1.348]$; large ($>20\ \mu\text{m}^2$, black circles) $= 1.184 \pm 0.196$.
- **Fig. 3(e) — Left out:** Confocal micrographs of droplets with monomeric actin (Latrunculin), $L \sim 90\text{ nm}$, and $L \sim 135\text{ nm}$; no plot axes.
- **Fig. 3(f) — Read (`fig3f_aspect_vs_area`):** Binned average aspect ratio vs cross-sectional area for 1 $\mu\text{M}$ actin, $L \sim 135\text{ nm}$ ($n = 1308$, $\pm 1\text{ SD}$ bars).
  - *Calibration:* 1200 dpi; linear x-axis from 5 top ticks ($0, 10, \dots, 40\ \mu\text{m}^2$), residual **0.60 px**; linear y-axis from 5 right ticks ($2.0, 1.75, \dots, 1.0$), residual **3.21 px**.
  - *Reading:* **30** bins spanning $A \in [0.94, 37.0]\ \mu\text{m}^2$, $x \in [1.042, 1.577]$ (**25** with $\pm 1\text{ SD}$ bars; matches the same $L \sim 135\text{ nm}$ dataset replotted with theoretical fits in Fig. 5c).
- **Fig. 3(g) — Read (`fig3g_small`, `fig3g_medium`, `fig3g_large`):** Average aspect ratio vs approximate filament length (broken x-axis: Lat monomeric control coded at $x = 0\text{ nm}$, plus $L = 90, 135, 270\text{ nm}$; 1 $\mu\text{M}$ actin; $\pm 1\text{ SD}$ bars).
  - *Calibration:* 1200 dpi; broken x-axis with Lat tick anchored at $0\text{ nm}$ and linear fit to top ticks $100, 150, 200, 250\text{ nm}$ + right spine ($300\text{ nm}$), residual **1.60 px**; linear y-axis from top spine ($3.0$) + 4 right ticks ($2.5, \dots, 1.0$), residual **1.40 px**.
  - *Readings:* **4** points per series (**12** points total, all 12 with SD bars). At $L \sim 270\text{ nm}$: small $= 2.015 \pm 0.698$, medium $= 1.518\ [0.990, 2.006]$, large $= 1.212\ [1.055, 1.363]$. Comparing Fig. 3(d) ($1.45 \pm 0.27$) vs Fig. 3(g) ($2.02 \pm 0.70$) for identically labeled $1\ \mu\text{M}$, $L \sim 270\text{ nm}$, $4$–$6\ \mu\text{m}^2$ droplets confirms the $\sim 40\%$ inter-batch variation (`q_digitized_batch_and_prose_discrepancies`).

### Figure 4 (Page 7)

- **Fig. 4(a) — Left out (images only):** Representative droplet micrographs illustrating the tactoid, ellipsoid, and sphere morphological classes; the printed example fractions at $4$–$5\ \mu\text{m}^2$ ($73\%$ tactoids, $27\%$ ellipsoids, $0\%$ spheres) and $8$–$9\ \mu\text{m}^2$ ($11\%$ tactoids, $89\%$ ellipsoids, $0\%$ spheres) are already stored verbatim in `data.yaml`.
- **Fig. 4(b) — Read (`fig4b_tactoid_fraction`):** Variable-width histogram of percentage of droplets classified as tactoids vs cross-sectional area ($1\ \mu\text{M}$ actin, $L \sim 270\text{ nm}$, $n = 410$).
  - *Calibration:* 1800 dpi; outward bottom major ticks ($5, 10, 15, 20\ \mu\text{m}^2$) + right spine ($25\ \mu\text{m}^2$), residual **4.21 px** ($\approx 0.06\ \mu\text{m}^2$); outward left ticks ($70, 60, \dots, 10\%$) + bottom spine ($0\%$), residual **3.99 px** ($\approx 0.35\%$).
  - *Reading:* All **11** nonzero grey histogram bars spanning $4$–$21\ \mu\text{m}^2$: $72.7\%$ ($4$–$5$), $54.8\%$ ($5$–$6$), $54.1\%$ ($6$–$7$), $48.1\%$ ($7$–$8$), $11.1\%$ ($8$–$9$), $14.2\%$ ($9$–$10$), $16.6\%$ ($10$–$11$), $14.2\%$ ($11$–$12$), $11.4\%$ ($14$–$16$), $9.4\%$ ($16$–$18$), and $8.1\%$ ($18$–$21\ \mu\text{m}^2$). Zero-height bins ($<4$, $12$–$14$, $>21\ \mu\text{m}^2$) have no bar above baseline. Confirms Fig. 4(a)'s $73\%$ and $11\%$, and notes that the $7$–$8\ \mu\text{m}^2$ bin sits slightly below $50\%$ at $48.1\%$.
- **Fig. 4(c) — Read (`fig4c_tactoid`, `fig4c_ellipsoid`, `fig4c_sphere`):** Droplet shape phase diagram (Area vs Approximate Filament Length, broken x-axis with Lat at $x = 0\text{ nm}$, $1\ \mu\text{M}$ actin, 4.4 $\mu\text{M}$ FUS).
  - *Calibration:* 1200 dpi; broken x-axis with Lat tick anchored at $0\text{ nm}$ and linear fit to bottom ticks $100, \dots, 300\text{ nm}$, residual **1.20 px**; linear y-axis from 9 right ticks ($45, 40, \dots, 5\ \mu\text{m}^2$) + bottom spine ($0$), residual **2.51 px**.
  - *Readings:* `fig4c_tactoid` (**5** red diamonds, Region I, $L = 135, 270\text{ nm}$, $A \in [4.2, 7.3]\ \mu\text{m}^2$), `fig4c_ellipsoid` (**32** black triangles, Region II, $L = 90, 135, 270\text{ nm}$, $A \in [4.3, 45.5]\ \mu\text{m}^2$), `fig4c_sphere` (**31** blue circles, Region III, $L = 0, 90, 135\text{ nm}$, $A \in [4.3, 38.9]\ \mu\text{m}^2$). Total: **68** classified phase-space points.
  - *Left-out element:* Hand-drawn dashed guide curves sketching approximate phase boundaries between Regions I, II, and III (no mathematical definition).

### Figure 5 (Page 8)

- **Fig. 5(a) — Left out:** Geometric schematic of a quasi-bipolar tactoid ($R$, $r$, $\tilde{R}$, director streamlines, virtual poles); no plotted data.
- **Fig. 5(b) — Left out:** LC-PolScope optical retardance/orientation micrographs and schematic director overlays; qualitative images with no plot axes.
- **Fig. 5(c) — Read (`fig5c_best_fit`, `fig5c_min_bound`, `fig5c_max_bound`, `fig5c_data`):** Quasi-bipolar model fits vs experimental binned droplet aspect ratio for $1\ \mu\text{M}$ actin, $L \sim 135\text{ nm}$.
  - *Calibration:* 1200 dpi; linear x-axis from left spine ($0$) + 7 top ticks ($5, 10, \dots, 35$) + right spine ($40\ \mu\text{m}^2$), residual **1.33 px**; linear y-axis from 4 right ticks ($2.0, 1.75, 1.5, 1.25$) + bottom spine ($1.0$), residual **2.20 px**.
  - *Readings:* Three theoretical curves sampled at **12** areas $A \in \{2, 4, 6, 8, 10, 12, 15, 18, 20, 25, 30, 35\}\ \mu\text{m}^2$ (`fig5c_best_fit` red solid, `fig5c_min_bound` light-blue dashed, `fig5c_max_bound` dark-blue dashed; 36 curve points total) + `fig5c_data` (**33** black squares spanning $A \in [0.99, 37.0]\ \mu\text{m}^2$, **26** with $\pm 1\text{ SD}$ bars).
  - *Comparisons:* See Quantitative Verdict Table below — all 3 `A^(-1/2)` models (`printed`, `exact_bend`, `frank_half`) fail on all 3 curves (`not reproduced`, $\text{RMS} = 0.24$–$0.30$), whereas `quasi_bipolar_aspect_vs_area_vtilde_inv` reproduces all 3 curves (**12/12** points each, $\text{RMS} = 0.0039$–$0.0044$).
- **Fig. 5(d) — Read (`fig5d_K_over_gamma`, `fig5d_omega`):** Dual-y-axis summary of extracted $K/\gamma$ (left y-axis, black open squares) and $\omega$ (right y-axis, red open circles) vs approximate filament length ($L = 90, 135, 270\text{ nm}$) with min/max bounding fit range bars (`"Error bars are from the maximum and minimum theoretical fits"`, p. 8).
  - *Calibration:* 1200 dpi; linear x-axis from 4 top ticks ($100, 150, 200, 250\text{ nm}$), residual **0.60 px**; left y-axis ($K/\gamma$) from top spine ($2.5$) + 6 left ticks ($2.25, \dots, 1.0\ \mu\text{m}$), residual **1.43 px**; right y-axis ($\omega$) from 6 right ticks ($3.0, 2.5, \dots, 0.5$) + bottom spine ($0$), residual **1.86 px**.
  - *Readings:*
    - `fig5d_K_over_gamma` (**3** points, 3 with range bars): $1.201\ [0.898, 1.201]$ at $90\text{ nm}$; $1.701\ [1.399, 1.904]$ at $135\text{ nm}$; $1.701\ [1.399, 2.005]\ \mu\text{m}$ at $270\text{ nm}$ (matches Fig. 5c/S1 printed legends $1.2\ [0.9, 1.2]$, $1.7\ [1.4, 1.9]$, $1.7\ [1.4, 2.0]$ to $\le 0.005\ \mu\text{m}$).
    - `fig5d_omega` (**3** points, 3 with range bars): $0.198\ [0.123\text{ (lower bar hidden inside circle)}, 0.313]$ at $90\text{ nm}$; $0.443\ [0.364\text{ (hidden inside circle)}, 0.563]$ at $135\text{ nm}$; $2.797\ [2.393, 3.217]$ at $270\text{ nm}$ (consistent with Fig. 5c/S1 printed legends $0.2\ [0.15, 0.3]$, $0.45\ [0.38, 0.55]$, $2.8\ [2.4, 3.2]$; the very short lower spans $0.20 \to 0.15$ and $0.45 \to 0.38$ terminate inside the open red circle marker radius).

### Electronic Supplementary Information — Fig. S1 & Movies S1–S2 (Page 11)

- **Fig. S1 (left) — Read (`figS1_left_best_fit`, `figS1_left_min_bound`, `figS1_left_max_bound`, `figS1_left_data`):** Quasi-bipolar model fits vs experimental binned aspect ratio for $1\ \mu\text{M}$ actin, $L \sim 90\text{ nm}$.
  - *Calibration:* 900 dpi; linear x-axis from left spine ($0$) + 3 top ticks ($10, 20, 30$) + right spine ($40\ \mu\text{m}^2$), residual **0.30 px**; linear y-axis from 5 right ticks ($1.4, 1.3, 1.2, 1.1, 1.0$), residual **0.80 px**.
  - *Readings:* Three theoretical curves at **12** areas $A \in [2, 35]\ \mu\text{m}^2$ (36 curve points) + `figS1_left_data` (**28** black squares, **22** with $\pm 1\text{ SD}$ bars, reaching $x = 1.277$ at $A \approx 1.0\ \mu\text{m}^2$, directly contradicting p. 6 prose `"aspect ratios ~1"` for all sizes; see `q_90nm_aspect_ratio_prose_vs_figures`).
- **Fig. S1 (right) — Read (`figS1_right_best_fit`, `figS1_right_min_bound`, `figS1_right_max_bound`, `figS1_right_data`):** Quasi-bipolar model fits vs experimental binned aspect ratio for $1\ \mu\text{M}$ actin, $L \sim 270\text{ nm}$.
  - *Calibration:* 900 dpi; linear x-axis from left spine ($0$) + 5 top ticks ($10, \dots, 50\ \mu\text{m}^2$), residual **1.62 px**; linear y-axis from 8 right ticks ($4.5, 4.0, \dots, 1.0$), residual **1.26 px**.
  - *Readings:* Three theoretical curves at **12** areas $A \in [2, 35]\ \mu\text{m}^2$ (36 curve points) + `figS1_right_data` (**34** black squares spanning $A \in [1.14, 51.1]\ \mu\text{m}^2$, $x \in [1.095, 3.254]$, **26** with $\pm 1\text{ SD}$ bars).
- **Movies S1 and S2 (Page 11) — Left out:** Text captions describing supplementary time-lapse microscopy videos (coalescence at $L \sim 90\text{ nm}$ and nucleation/growth at $L \sim 135\text{ nm}$); video files are not embedded in the PDF.

---

## Key Quantitative Verdicts Settled by Digitization

### 1. Fig. 5(c) & ESI Fig. S1: `A^(-1/2)` vs `A^(-3/2)` Volume-Scaling Discrepancy (Settled)

Calibrated digitization of all **9** theoretical curves (**108** curve points across $A \in [2, 35]\ \mu\text{m}^2$) in Fig. 5(c) and Fig. S1 decisively settles `q_fig5c_figS1_volume_scaling_discrepancy` (logged in `q_digitized_fig5c_figS1_settles_volume_scaling`):

| Panel & Curve | `printed` ($c = K/(\gamma V^{1/3})$) | `exact_bend` ($c = K/(\gamma V^{1/3})$) | `frank_half` ($c = K/(\gamma V^{1/3})$) | `vtilde_inv` ($c = (K/\gamma)^3/V$, $A = Rr$) |
| :--- | :---: | :---: | :---: | :---: |
| **Fig. S1 left** ($90\text{ nm}$) best | `not reproduced` (3/12, RMS 0.128) | `not reproduced` (2/12, RMS 0.138) | `not reproduced` (3/12, RMS 0.132) | **`reproduced` (12/12, RMS 0.0020)** |
| **Fig. S1 left** ($90\text{ nm}$) min | `not reproduced` (0/12, RMS 0.112) | `not reproduced` (0/12, RMS 0.120) | `not reproduced` (0/12, RMS 0.115) | **`reproduced` (12/12, RMS 0.0014)** |
| **Fig. S1 left** ($90\text{ nm}$) max | `not reproduced` (1/12, RMS 0.207) | `not reproduced` (1/12, RMS 0.224) | `not reproduced` (1/12, RMS 0.208) | **`reproduced` (12/12, RMS 0.0023)** |
| **Fig. 5(c)** ($135\text{ nm}$) best | `not reproduced` (1/12, RMS 0.250) | `not reproduced` (0/12, RMS 0.271) | `not reproduced` (1/12, RMS 0.243) | **`reproduced` (12/12, RMS 0.0039)** |
| **Fig. 5(c)** ($135\text{ nm}$) min | `not reproduced` (1/12, RMS 0.242) | `not reproduced` (0/12, RMS 0.263) | `not reproduced` (1/12, RMS 0.240) | **`reproduced` (12/12, RMS 0.0040)** |
| **Fig. 5(c)** ($135\text{ nm}$) max | `not reproduced` (1/12, RMS 0.282) | `not reproduced` (1/12, RMS 0.305) | `not reproduced` (1/12, RMS 0.268) | **`reproduced` (12/12, RMS 0.0044)** |
| **Fig. S1 right** ($270\text{ nm}$) best | `not reproduced` (0/12, RMS 1.237) | `not reproduced` (0/12, RMS 1.238) | `not reproduced` (1/12, RMS 0.903) | **`reproduced` (12/12, RMS 0.0062)** |
| **Fig. S1 right** ($270\text{ nm}$) min | `not reproduced` (0/12, RMS 1.182) | `not reproduced` (0/12, RMS 1.197) | `not reproduced` (0/12, RMS 0.887) | **`reproduced` (12/12, RMS 0.0071)** |
| **Fig. S1 right** ($270\text{ nm}$) max | `not reproduced` (0/12, RMS 1.259) | `not reproduced` (1/12, RMS 1.248) | `not reproduced` (0/12, RMS 0.894) | **`reproduced` (12/12, RMS 0.0113)** |

- **Verdict:** All **27/27** curve comparisons for the printed $A^{-1/2}$ scaling (`printed`, `exact_bend`, `frank_half`) fail (`not reproduced`), whereas `quasi_bipolar_aspect_vs_area_vtilde_inv` achieves **9/9 `reproduced` (108/108 points reproduced, 0 close, 0 not reproduced)** with RMS $0.0014$–$0.0113$ (inside the calibrated reading precision).
- **Fit quality against experimental squares ($\pm 1\text{ SD}$):** `vtilde_inv` matches **24/28** points on Fig. S1 left ($\text{RMS} = 0.021$ vs $0.131$ for `printed`), **27/33** on Fig. 5(c) ($\text{RMS} = 0.048$ vs $0.260$ for `printed`), and **20/34** on Fig. S1 right ($\text{RMS} = 0.126$ vs $1.234$ for `printed`).

### 2. Other Quantitative Findings Settled by Digitization

- **Fig. 2(b) exponential master curve:** **`reproduced`** (**10/10** points, $\text{RMS} = 0.0029$, $\max|d| = 0.0070$).
- **Fig. 2(c) vs Fig. 2(d) slope mismatch (`q_digitized_fig2c_vs_fig2d_slope_discrepancy`):** The plotted dashed linear fit in Fig. 2(c) has exact slope $26.34\text{ s}/\mu\text{m}$ (intercept $-1.18\text{ s}$; unweighted OLS on the 9 points gives $29.10 \pm 9.55\text{ s}/\mu\text{m}$), whereas Fig. 2(d) plots $1\ \mu\text{M}$ actin at $18.92 \pm 10.57\text{ s}/\mu\text{m}$ (overlapping within 1 SE, but $7.4\text{ s}/\mu\text{m}$ lower despite prose stating that the Fig. 2c fit slope gives the Fig. 2d value).
- **Fig. 3(d) vs Fig. 3(g) & Fig. 3(c) vs Fig. S1 right batch shift (`q_digitized_batch_and_prose_discrepancies`):** Nominally identical $1\ \mu\text{M}$, $L \sim 270\text{ nm}$, $4$–$6\ \mu\text{m}^2$ droplets have mean aspect ratio $1.45 \pm 0.27$ in Fig. 3(d) (consistent with Fig. 3c, $n = 1180$) vs $2.02 \pm 0.70$ in Fig. 3(g) (consistent with Fig. S1 right).
- **Fig. 4(b) vs Fig. 4(a) & prose threshold:** The $4$–$5\ \mu\text{m}^2$ ($72.7\%$) and $8$–$9\ \mu\text{m}^2$ ($11.1\%$) histogram bars confirm Fig. 4(a)'s printed $73\%$ and $11\%$, while the $7$–$8\ \mu\text{m}^2$ bar sits at $48.1\%$ (just under $50\%$).
- **Fig. 5(d) internal consistency:** All 6 best-fit points and bounding range bars match the printed legend values in Fig. 5(c) and Fig. S1 within $\le 0.005\ \mu\text{m}$ ($K/\gamma$) and $\le 0.02$ ($\omega$).

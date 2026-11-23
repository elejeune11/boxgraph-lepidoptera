# Figure Inventory and Digitization Log — `sekine2024`

Complete panel-by-panel inventory of `corpus/actin_bundling/2024_Sekine.pdf` (38 PDF pages: pp. 1–13 main article & Methods, pp. 14–27 Supplementary Figs. 1–7 & Supplementary Tables 1–6, pp. 28–36 Supplementary Notes & Supplementary References, pp. 37–38 Description of Additional Supplementary Files). Note that the published article contains **Main Figs. 1–5** (there is no Fig. 6), **Supplementary Figs. 1–7** (there is no Supplementary Fig. 8), and **SI Notes Fig. S1** (p. 29).

All **44 quantitative readings** in [`figures.yaml`](figures.yaml) were implemented independently from the PDF as 900–1200 dpi `mode: code` Python scripts under [`figures/`](figures/), verified with `boxgraph figure check`, and visually inspected via every generated magenta overlay PNG in `figwork/`. No files from `sekine_source_data/` were inspected.

---

## Page 1 (printed p. 1)
- No figures (Title, Abstract, Introduction text only).

## Page 2 (printed p. 2) — Main Figure 1: Emergence of actin nanoclusters during tracheal tube expansion
- **Fig. 1a**: Left out — schematic of the tracheal tube lumen/membrane paired with a representative phalloidin F-actin fluorescence micrograph (scale bar 5 µm; no quantitative plot axes).
- **Fig. 1b**: Left out — qualitative developmental timeline schematic (`-80` to `20 min`, Stage 14–16) with schematic gradient bars and no quantitative y-axis.
- **Fig. 1c**: Left out — confocal fluorescence micrographs and binarized inset crops of fixed Stage 15 `btl>lifeact::GFP` tracheal tubes.
- **Fig. 1d**: Left out — confocal fluorescence micrograph and binarized inset crop of fixed Early Stage 16 `btl>lifeact::GFP` tracheal tube.
- **Fig. 1e**: Left out — polar probability histogram of cluster major-axis orientation (`0°`–`90°`, St 15 vs St 16) whose radial `Probability` axis carries concentric grid rings with **no numerical radial tick labels**, making the radial scale uncalibratable from the figure alone (circular summary medians/means are already recorded from Supplementary Table 1 in `data.yaml`).
- **Fig. 1f (Aspect ratio density, Stage 15 vs Stage 16)**: Digitized both smooth kernel density curves at 1200 dpi (`rx = 1.04 px`, `ry = 0.34 px`, 12 points each at `AR = 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0`):
  - `fig1f_st15` (`d_fig1f_st15`, red St 15 curve, peak density `0.806` at `AR = 1.4`, decaying to `~0.00` by `AR = 6.0`).
  - `fig1f_st16` (`d_fig1f_st16`, teal St 16 curve, broader shifted distribution peaking at `0.522` at `AR = 1.6` with elevated tail `0.019` at `AR = 6.0`).
- **Fig. 1g**: Left out — live time-lapse confocal fluorescence micrographs (`-60` to `20 min`, scale bar 2 µm).
- **Fig. 1h (Actin orientation vs time)**: Digitized as `fig1h_orientation` (`d_fig1h_orientation`, 1200 dpi, `rx = 0.39 px`, `ry = 0.36 px`, 11 points `-80` to `20 min`, `error_bars: sd`). Mean dominant direction rises from `15.65°` (`-80 min`) to `83.90°` (`0*`–`20 min`) with SD ribbon `[3.48°..86.95°]`.
- **Fig. 1i (Actin connectivity / coherency vs time)**: Digitized as `fig1i_coherency` (`d_fig1i_coherency`, 1200 dpi, `rx = 0.97 px`, `ry = 0.60 px`, 11 points `-80` to `20 min`, `error_bars: sd`). Mean coherency rises from `0.1345` (`-80 min`) to `0.4762` (`20 min`) with SD ribbon `[0.0884..0.5588]`. *(Minor printed y-axis typo noted: `"Ccoherency (0-1)"`.)*
- **Fig. 1j (Tube expansion / diameter vs time)**: Digitized as `fig1j_diameter` (`d_fig1j_diameter`, 1200 dpi, `rx = 0.69 px`, `ry = 0.72 px`, 6 points `-80` to `20 min` at 20-min steps, `error_bars: range`). Mean dorsal-trunk diameter rises from `2.598 µm` to `3.942 µm` (min–max ribbon `[2.219..4.655 µm]`). *(Minor printed y-axis typo noted: `"Tube diamter (µm)"`.)*
- **Fig. 1k (Tube elongation / length vs time)**: Digitized as `fig1k_length` (`d_fig1k_length`, 1200 dpi, `rx = 0.54 px`, `ry = 0.44 px`, 6 points `-80` to `20 min`, `error_bars: range`). Mean dorsal-trunk length rises from `16.73 µm` to `24.24 µm` (min–max ribbon `[13.55..26.91 µm]`).

## Page 3 (printed p. 3)
- No figures (Results text only).

## Page 4 (printed p. 4) — Main Figure 2: Circumferential stabilization of nanocluster fusion
- **Fig. 2a**: Left out — live fluorescent micrographs, binarized fusion-site masks, and spectrum-coded temporal projections (`-40 min` and `0* min`, scale bar 0.5 µm).
- **Fig. 2b (Nanocluster fusion duration vs time and fusion angle)**: Digitized the 15 plotted black-diamond group means across the 5 developmental time points (`-40, -30, -20, -10, 0* min`) at 900 dpi (`ry = 0.61 px`). Individual jittered discrete-integer frame counts (`1..11` frames, `n = 74..261` per box) heavily overlap at identical integer rows and share discrete median bars (`2` or `3` frames), so the continuous black-diamond means carry the primary quantitative signal:
  - `fig2b_mean_60_90` (`d_fig2b_mean_60_90`, circumferential `60–90°`, blue boxes): `2.88, 3.20, 3.31, 3.22, 3.57` frames (`-40` to `0* min`).
  - `fig2b_mean_30_60` (`d_fig2b_mean_30_60`, diagonal `30–60°`, pink boxes): `2.46, 2.69, 2.49, 2.81, 2.63` frames.
  - `fig2b_mean_0_30` (`d_fig2b_mean_0_30`, longitudinal `0–30°`, light-pink boxes): `2.50, 2.38, 2.24, 2.33, 2.62` frames.
- **Fig. 2c**: Left out — qualitative schematic cartoon of repetitive nanocluster fusion/fission and circumferential stabilization.

## Page 5 (caption on Page 6) — Main Figure 3: RNAi/mutant screen, colocalization, and fusion stability
- **Fig. 3a (RNAi/mutant screen: average major length vs median orientation at St. 16)**: Digitized all 16 highlighted colored markers at 900 dpi (`rx = 0.50 px`, `ry = 1.30 px`). Unlabeled mild-defect gray background circles are omitted here because all 28 conditions are already tabulated in Supplementary Table 2 (`data.yaml`) and digitized individually in Supplementary Fig. 3b:
  - `fig3a_controls` (`d_fig3a_controls`, 2 red control markers): `(63.49°, 363.4 nm)` (mislabeled `Control (attP40)` on the plot, matching Supplementary Table 2's `attP2`: `63.5°, 363.37 nm`) and `(64.52°, 345.6 nm)` (mislabeled `Control (attP2)`, matching Table 2's `attP40`: `64.5°, 345.57 nm`) — directly confirming `q_control_labels_swapped_fig3a`.
  - `fig3a_strong_rnai` (`d_fig3a_strong_rnai`, 10 blue strong-defect RNAi markers spanning `44.49°..60.49°` and `279.6..314.2 nm`).
  - `fig3a_mutants` (`d_fig3a_mutants`, 4 green mutant/DN markers spanning `42.54°..62.50°` and `261.4..360.0 nm`; notably `myosin II^DN` sits at `(54.42°, 325.6 nm)` vs `RhoA^DN` at `(54.53°, 279.0 nm)`, confirming Table 2's `54.4°` is real rather than a `.5` typo).
- **Fig. 3b**: Left out — confocal projection and binarized micrographs of early stage 16 control, RNAi, and mutant tracheal tubes.
- **Fig. 3c**: Left out — qualitative protein–protein interaction network diagram grouping screen-positive genes by function.
- **Fig. 3d (Stage 15 colocalization with actin, Pearson's $R$)**: Digitized all 7 protein boxplots (median + min–max whisker/outlier range) at 900 dpi (`ry = 0.12 px`, `error_bars: range`) as `fig3d_coloc_st15` (`d_fig3d_coloc_st15`): `1=GFP` `-0.115 [-0.335..0.095]`, `2=DAAM` `0.500 [0.415..0.736]`, `3=α-Actinin` `0.270 [0.035..0.387]`, `4=Zasp52` `0.720 [-0.420..0.886]` (capturing the lower black outlier dot at `-0.420`), `5=Myosin II` `0.240 [-0.143..0.487]`, `6=ROCK` `0.335 [-0.106..0.856]`, `7=RhoGEF2` `-0.160 [-0.235..-0.106]`.
- **Fig. 3e (Early stage 16 fusion duration by genotype and fusion angle)**: Digitized all 12 black-diamond group means across `1=Control, 2=DAAM RNAi, 3=zasp52 RNAi, 4=myosin II^DN` at 900 dpi (`ry = 0.09 px`):
  - `fig3e_mean_60_90` (`d_fig3e_mean_60_90`, `60–90°`): `3.41, 2.74, 3.50, 2.69` frames.
  - `fig3e_mean_30_60` (`d_fig3e_mean_30_60`, `30–60°`): `2.67, 2.26, 2.42, 2.08` frames.
  - `fig3e_mean_0_30` (`d_fig3e_mean_0_30`, `0–30°`): `2.57, 2.26, 2.18, 2.14` frames.
- **Fig. 3f**: Left out — confocal projection and binarized micrographs of late stage 16 control, RNAi, and mutant tracheal tubes.

## Page 7 (caption on Page 8) — Main Figure 4: Agent-based CGMD simulation results
- **Fig. 4a**: Left out — schematic diagram of coarse-grained actin filament, α-Actinin crosslinker, and Myosin II motor primitives.
- **Fig. 4b**: Left out — 4 × 3 grid (`1-L` through `4-H`) of CGMD simulation snapshots, continuous density maps, and binarized maps (qualitative spatial patterns; quantitative summary metrics are plotted in Fig. 4c and Supplementary Fig. 6b and tabulated in Supplementary Table 3).
- **Fig. 4c (CGMD + experimental average major length vs median orientation)**: Digitized all 9 unoccluded markers at 900 dpi (`rx = 0.91 px`, `ry = 2.74 px`) as `fig4c_sim_exp` (`d_fig4c_sim_exp`): `3-H` `(38.98°, 347.1 nm)`, `1-M` `(41.96°, 318.7 nm)`, `St.15` `(46.98°, 328.8 nm)`, `2-M` `(53.92°, 283.0 nm)`, `4-H` `(57.01°, 389.6 nm)`, `3-M` `(58.92°, 349.5 nm)`, `2-H` `(58.92°, 439.9 nm)`, `St.16` `(66.50°, 383.6 nm)`, and `4-M` `(66.92°, 412.0 nm)` — all matching Supplementary Table 3 within `0.1–0.6 nm`. Condition `1-H` (`39°, 344.5 nm`) is completely occluded directly behind `3-H` in Fig. 4c (and is cleanly resolved in Supplementary Fig. 6b/6c below).
- **Fig. 4d**: Left out — three dense overlapping hollow-circle scatter clouds of `(Δx, Δy)` single-step displacements (`n = 1461` in `1-M`, `n = 1134` in `2-M`, `n = 152` in `4-M`) that cannot be separated individually; their printed summary anisotropy ratios `|Δy|/|Δx|` (`0.98, 1.13, 1.45`) and component means/SDs are already recorded from Fig. 4d & Supplementary Table 5 and verified via `cmp_aniso_sim_*`.
- **Fig. 4e (CGMD shift magnitude boxplots `|Δx|` vs `|Δy|`)**: Digitized all 6 boxplots (median + whisker range, in µm) across the 3 subpanels (`1=1-M Δx, 2=1-M Δy, 3=2-M Δx, 4=2-M Δy, 5=4-M Δx, 6=4-M Δy`) at 1200 dpi (`ry = 1.31, 0.49, 0.76 px`, `error_bars: range`) as `fig4e_shift` (`d_fig4e_shift`): medians `0.0141, 0.0123, 0.0133, 0.0140, 0.0141, 0.0158 µm` with upper whiskers rising from `0.0439–0.0442 µm` (`1-M`) to `0.0571–0.0589 µm` (`4-M`) and lower whiskers reaching `~0.000–0.001 µm`.
- **Fig. 4f (4-M cable self-organization time course, `t = 0..40 s`)**: Digitized both 41-point time-series curves (Excel category-axis chart with 41 equal 1-s bins, `ry = 1.60 px`, 1200 dpi; bottom snapshots at `t = 0, 4, 40 s` are non-quantitative density renderings):
  - `fig4f_orientation` (`d_fig4f_orientation`): dominant direction jumps rapidly from initial random orientations (`4.66°`–`89.6°` at `t = 0–5 s`) and stabilizes near `77°–90°` (`89.30°` at `t = 40 s`).
  - `fig4f_coherency` (`d_fig4f_coherency`): starts at `0.0215` (`t = 0 s`), peaks at `0.2238` (`t = 31 s`), and finishes at `0.1891` (`t = 40 s`) — never approaching Supplementary Table 3's reported `0.30` (the `> 0.3` cable-formation criterion) within the plotted 40-s window, quantitatively corroborating item (3) of `q_fig5f_and_figs5d_significance_claims`.

## Page 9 (caption on Page 10) — Main Figure 5: *In vivo* nanocluster motion anisotropy & PIV
- **Fig. 5a**: Left out — five dense overlapping scatter clouds of `(Δx, Δy)` nanocluster displacements (`-40` to `0* min`, `n = 314..372` points per panel); their printed summary anisotropy ratios `|Δy|/|Δx|` (`1.08, 1.36, 1.31, 1.39, 1.72`) and component means/SDs are already recorded and verified from Supplementary Table 5 (`cmp_aniso_st15_*`).
- **Fig. 5b**: Left out — 2D color-coded PIV vector-field maps overlaid on live cells (`-40` to `0* min`, scale bar 0.5 µm); quantitative velocity averages are plotted in Fig. 5c.
- **Fig. 5c (Stage 15 PIV mean velocity `|Vx|` and `|Vy|` vs time)**: Digitized both mean ± SD time-course curves across `-40, -30, -20, -10, 0* min` at 1200 dpi (`rx = 0.69 px`, `error_bars: sd`):
  - `fig5c_vx` (`d_fig5c_vx`, coral-red longitudinal `|Vx|`): `27.5 ± 5.1, 20.8 ± 3.2, 22.2 ± 4.8, 21.3 ± 3.9, 20.7 ± 3.8 nm/s`.
  - `fig5c_vy` (`d_fig5c_vy`, teal circumferential `|Vy|`): `30.3 ± 4.9, 25.6 ± 2.9, 32.1 ± 4.8, 29.4 ± 2.1, 34.8 ± 0.8 nm/s`.
- **Fig. 5d**: Left out — four dense overlapping scatter clouds of `(Δx, Δy)` displacements at early stage 16 (`n = 352..977` points per genotype); their printed summary ratios `|Δy|/|Δx|` (`1.60, 0.96, 0.98, 1.20`) and means/SDs are already recorded and verified in `data.yaml` / `comparisons.yaml` (`cmp_aniso_st16_*`).
- **Fig. 5e**: Left out — 2D color-coded PIV vector-field maps across the four genotypes (`Control, DAAM RNAi, zasp52 RNAi, myosin II^DN`).
- **Fig. 5f (Early stage 16 PIV mean velocity boxplots `|Vx|` vs `|Vy|`)**: Digitized all 8 paired boxplots (median + full min–max dot range across `1=Control, 2=DAAM RNAi, 3=Zasp52 RNAi, 4=myosin II^DN`) at 1200 dpi (`ry = 1.11 px`, `error_bars: range`):
  - `fig5f_vx` (`d_fig5f_vx`, coral-pink `|Vx|`): `17.1 [10.2..24.3]`, `29.4 [18.3..33.3]`, `24.3 [19.4..38.6]`, `22.6 [16.9..26.2] nm/s`.
  - `fig5f_vy` (`d_fig5f_vy`, cyan `|Vy|`): `28.8 [27.9..37.7]`, `28.4 [24.4..32.7]`, `26.5 [24.5..28.8]`, `27.2 [22.4..29.5] nm/s`.
- **Fig. 5g**: Left out — conceptual three-circle Venn diagram classifying RNAi phenotypes across nanocluster self-organization, myosin II fusion, and anisotropy sensing.

## Pages 11–14 (printed pp. 11–13 + SI Cover p. 14)
- No figures (Methods, References, Author Information, and Supplementary Information title page). **No Main Fig. 6 exists in the paper.**

## Page 15 (caption continued on Page 16) — Supplementary Figure 1: Cluster sizes, cable spacing, RhoGEF2, and skeleton elongation
- **Supp. Fig. 1a**: Left out — confocal fluorescence and binarized micrographs of Stage 14 and Stage 17 tracheal tubes.
- **Supp. Fig. 1b (Major axis length density, St 15 vs St 16)**: Digitized both smooth KDE curves at 1200 dpi (`rx = 0.07 px`, `ry = 0.00 px`, 11 points each `0.1..1.0 µm`):
  - `figs1b_st15` (`d_figs1b_st15`, red St 15 curve, peak `3.246 µm^-1` at `0.20 µm`).
  - `figs1b_st16` (`d_figs1b_st16`, teal St 16 curve, broader peak `2.598 µm^-1` at `0.25 µm` with longer tail `0.153 µm^-1` at `1.0 µm`).
- **Supp. Fig. 1c (Minor axis length density, St 15 vs St 16)**: Digitized both smooth KDE curves at 1200 dpi (`rx = 0.06 px`, `ry = 0.06 px`):
  - `figs1c_st15` (`d_figs1c_st15`, red St 15 curve, 11 unoccluded points across `0.10..0.60 µm`, peaking at `5.896 µm^-1` at `0.16 µm`; points around `0.08, 0.14–0.15, 0.18, 0.30, 0.45 µm` are completely covered by the overlapping teal St 16 stroke drawn on top).
  - `figs1c_st16` (`d_figs1c_st16`, teal St 16 curve, 11 points across `0.08..0.45 µm`, peaking at `6.217 µm^-1` at `0.12 µm`).
- **Supp. Fig. 1d (WT early stage 16 periodic cable interval boxplot)**: Digitized as `figs1d_intervals` (`d_figs1d_intervals`, 1200 dpi, `ry = 0.17 px`, `error_bars: range`): median `0.4468 µm` (`446.8 nm`), min–max range `[0.4352..0.4746 µm]` across `N = 8` embryos (`n = 294` cables), closely matching the main-text mean `450 ± 13.2 nm`.
- **Supp. Fig. 1e**: Left out — fluorescence micrographs of Phalloidin555, Utrophin::GFP, and lifeact::mScarletx2 at stages 15 and 16.
- **Supp. Fig. 1f (Junctional RFP::RhoGEF2 longitudinal/circumferential intensity ratio)**: Digitized both boxplots at 1200 dpi (`ry = 0.55 px`, `error_bars: range`) as `figs1f_rhogef2` (`d_figs1f_rhogef2`): `1=Stage 15 (n=5)` median `0.8447` (range `[0.760..1.212]`), `2=Stage 16 (n=5)` median `0.7692` (range `[0.701..1.022]`).
- **Supp. Fig. 1g (Actin elongation — mean skeleton length vs time)**: Digitized as `figs1g_skeleton` (`d_figs1g_skeleton`, 1200 dpi, `rx = 1.27 px`, `ry = 0.00 px`, 11 points `-80` to `20 min`, `error_bars: sd`). Mean skeleton length rises from `0.2923 µm` (`-80 min`) to `0.6510 µm` (`20 min`) with SD ribbon `[0.2263..0.9394 µm]`.
- **Supp. Fig. 1h**: Left out — transverse (`yz`) confocal cross-section micrographs of tracheal tube expansion (`-80` to `20 min`).

## Page 17 — Supplementary Figure 2: Fusion-site watershed detection example
- **Supp. Fig. 2a**: Left out — zoomed binarized micrograph showing watershed segmentation of candidate fusion sites (`#1`–`#4`).
- **Supp. Fig. 2b**: Left out — 6-frame binarized image montage illustrating one 4-frame (`1.28 s`) fusion event (printed example values already captured in `data.yaml` and verified via `cmp_fusion_dur_4frames`).

## Page 18 — Supplementary Figure 3: Quantitative RNAi screen summary
- **Supp. Fig. 3a**: Left out — 23 heavily overlapping thin OrientationJ distribution curves crossing densely throughout `0°–90°`, where individual gene traces cannot be unambiguously separated at pixel level (the median orientation of every single line is already recorded in Supplementary Table 2).
- **Supp. Fig. 3b (Average major length bar chart across 27 genetic conditions)**: Digitized all 27 blue bars in plotted left-to-right order at 900 dpi (`ry = 1.95 px`) as `figs3b_major_bars` (`d_figs3b_major_bars`, spanning `261.4 nm` for bar 25 `RhoA RNAi` to `382.0 nm` for bar 9 `twf`; omits `Myosin II[DN]` from Table 2's 28 rows). Crucially, Bar 1 `(Control) attP40 = 346.5 nm` is lower than Bar 2 `(Control) attP2 = 364.0 nm`, matching Supplementary Table 2 (`345.57` vs `363.37 nm`) and confirming that Fig. 3a and Supp. Fig. 3c swapped the red control labels.
- **Supp. Fig. 3c (Average major vs minor length scatter)**: Digitized the 2 highlighted red control markers at 900 dpi (`rx = 0.20 px`, `ry = 1.40 px`) as `figs3c_controls` (`d_figs3c_controls`): `(165.8 nm, 363.6 nm)` (mislabeled `Control (attP40)`, matching Table 2's `attP2`: `166 nm, 363.37 nm`) and `(160.9 nm, 348.7 nm)` (mislabeled `Control (attP2)`, matching Table 2's `attP40`: `161.02 nm, 345.57 nm`). Unlabeled light-blue RNAi markers are partially overlapping and are already tabulated individually in Supplementary Table 2 (`data.yaml`).

## Page 19 — Supplementary Figure 4: MIST v5.0 interaction network
- **Supp. Fig. 4**: Left out — database protein–protein and genetic interaction network graph (MIST v5.0); no numerical plot axes.

## Page 20 — Supplementary Figure 5: RhoA phenotypes & Stage 16 colocalization
- **Supp. Fig. 5a**: Left out — confocal projection and binarized micrographs of `RhoA RNAi` tracheal cells.
- **Supp. Fig. 5b**: Left out — confocal projection and binarized micrographs of `RhoA^DN` tracheal cells.
- **Supp. Fig. 5c**: Left out — two-channel confocal fluorescence micrographs of 6 tagged proteins vs actin at stages 15 and 16.
- **Supp. Fig. 5d (Stage 16 colocalization with actin, Pearson's $R$)**: Digitized all 7 protein boxplots (median + min–max whisker/outlier range) at 900 dpi (`ry = 0.11 px`, `error_bars: range`) as `figs5d_coloc_st16` (`d_figs5d_coloc_st16`): `1=GFP` `0.045 [-0.269..0.280]`, `2=DAAM` `0.030 [-0.366..0.246]`, `3=α-Actinin` `0.151 [-0.096..0.445]`, `4=Zasp52` `0.530 [0.463..0.727]`, `5=MyosinII` `0.250 [0.078..0.617]`, `6=ROCK` `0.190 [-0.285..0.556]`, `7=RhoGEF2` `-0.286 [-0.357..-0.049]` (capturing the upper black outlier dot at `-0.049`).

## Page 21 — Supplementary Figure 6: Molecular packing schematic & CGMD major vs minor lengths
- **Supp. Fig. 6a**: Left out — schematic cartoon of lifeact::GFP-labeled F-actin and α-Actinin/Zasp52 crosslinked dimer packing (its printed molecular dimensions `7, 17, 35, 58, 43, 100, 180 nm` are already recorded and checked via `cmp_width_single` and `cmp_width_dimer`).
- **Supp. Fig. 6b (Fig. 4b CGMD + experimental average major vs minor length)**: Digitized all 10 markers at 1200 dpi (`rx = 3.70 px`, `ry = 2.40 px`) as `figs6b_major_minor` (`d_figs6b_major_minor`), resolving all 8 CGMD conditions including `1-H` unoccluded at `(210.2 nm, 344.8 nm)` (vs Table 3: `210.1, 344.5 nm`), alongside `St.16` `(168.9, 383.5)`, `St.15` `(178.5, 328.7)`, `2-M` `(190.0, 282.1)`, `3-H` `(192.7, 347.5)`, `4-H` `(203.8, 388.9)`, `3-M` `(204.4, 351.1)`, `1-M` `(204.8, 317.9)`, `2-H` `(207.0, 439.1)`, and `4-M` `(220.7, 411.5 nm)`.
- **Supp. Fig. 6c (Supp. Fig. 7a CGMD average major vs minor length)**: Digitized all 7 markers (`1-M`, `1-H` circles + `a3-1, a4-1, a5-2, a6-1, a6-3` X-crosses) at 1200 dpi (`rx = 1.65 px`, `ry = 1.05 px`) as `figs6c_major_minor` (`d_figs6c_major_minor`): `(204.5, 318.3)`, `(206.5, 257.6)`, `(209.9, 345.2)`, `(213.1, 271.4)`, `(214.7, 317.5)`, `(218.8, 310.1)`, `(221.2, 317.9) nm`.

## Page 22 (caption continued on Page 23) — Supplementary Figure 7: Parameter sweeps & static structure factors $S(q)$
- **Supp. Fig. 7a (a1–a6)**: Left out — CGMD parameter-sweep simulation snapshots without anisotropy (`a1`–`a6`, 18 snapshots showing cluster morphology vs filament/crosslinker/motor parameters; quantitative sizes for selected conditions are digitized in Supp. Fig. 6c and tabulated in Supp. Table 3).
- **Supp. Fig. 7b (b1–b7)**: Left out — CGMD parameter-sweep simulation snapshots with friction anisotropy (`b1`–`b7`, 25 snapshots showing stripe/labyrinth/cluster transitions).
- **Supp. Fig. 7c1 (`1-M` clusters $S(q)$ vs $q/(2\pi)$)**: Digitized at 1200 dpi (`rx = 0.71 px`, `ry = 0.01 px`):
  - `figs7c1_fit` (`d_figs7c1_fit`, 9 points `q/(2π) = 1.5..3.1 µm^-1` on black dotted Gaussian fit): `cmp_sf_freq_c1_fit` (`gaussian_structure_factor_freq_variant`) **`reproduced`** (`med_abs = 0.0040`, `max_abs = 0.0110`, `max_rel = 1.81%`), whereas `cmp_sf_literal_c1_fit` (`gaussian_structure_factor_literal_p11`) **`not reproduced`** (`max_rel = 100.0%`).
  - `figs7c1_sim` (`d_figs7c1_sim`, all 30 blue simulated markers on uniform bin centers `q/(2π) = 0.1, 0.3, ..., 5.9 µm^-1`, peak `1.422` at `2.3 µm^-1`, high-$q$ baseline `~0.11`): `cmp_sf_freq_c1_sim` (`fit_quality`: `not reproduced` over the full 0.1–5.9 range because the Gaussian fits only the main peak and decays to zero outside it).
- **Supp. Fig. 7c2 (`3-M` labyrinth $S(q)$ vs $q/(2\pi)$)**: Digitized at 1200 dpi (`rx = 0.65 px`, `ry = 0.69 px`):
  - `figs7c2_fit` (`d_figs7c2_fit`, 6 points `1.1..2.1 µm^-1`): `cmp_sf_freq_c2_fit` **`reproduced`** (`med_abs = 0.0046`, `max_abs = 0.0129`, `max_rel = 2.84%`) vs `cmp_sf_literal_c2_fit` **`not reproduced`** (`min_rel = 29.5%`, `max_rel = 100.0%`).
  - `figs7c2_sim` (`d_figs7c2_sim`, 30 bins `0.1..5.9 µm^-1`, peak `0.773` at `1.9 µm^-1`): `cmp_sf_freq_c2_sim` (`fit_quality`: `not reproduced` over full spectrum).
- **Supp. Fig. 7c3 (`4-M` stripes $S(q)$ vs $q/(2\pi)$)**: Digitized at 1200 dpi (`rx = 0.67 px`, `ry = 0.46 px`):
  - `figs7c3_fit` (`d_figs7c3_fit`, 8 points `1.1..2.5 µm^-1`): `cmp_sf_freq_c3_fit` **`reproduced`** (`med_abs = 0.0015`, `max_abs = 0.0065`, `max_rel = 1.65%`) vs `cmp_sf_literal_c3_fit` **`not reproduced`** (`max_rel = 100.0%`).
  - `figs7c3_sim` (`d_figs7c3_sim`, 30 bins `0.1..5.9 µm^-1`, peak `0.473` at `1.7 µm^-1`): `cmp_sf_freq_c3_sim` (`fit_quality`: `not reproduced` over full spectrum).
- **Supp. Fig. 7c4 (`a4-1` low-turnover clusters $S(q)$ vs $q/(2\pi)$)**: Digitized at 1200 dpi (`rx = 0.57 px`, `ry = 0.64 px`):
  - `figs7c4_fit` (`d_figs7c4_fit`, 9 points `1.1..2.7 µm^-1`): `cmp_sf_freq_c4_fit` **`reproduced`** (`med_abs = 0.0091`, `max_abs = 0.0320`, `max_rel = 1.78%`) vs `cmp_sf_literal_c4_fit` **`not reproduced`** (`max_rel = 100.0%`).
  - `figs7c4_sim` (`d_figs7c4_sim`, 30 bins `0.1..5.9 µm^-1`, peak `3.068` at `2.1 µm^-1`): `cmp_sf_freq_c4_sim` (`fit_quality`: `not reproduced` over full spectrum).

**No Supplementary Fig. 8 exists in the PDF.**

## Pages 24–38 — Supplementary Tables 1–6, Supplementary Notes (Fig. S1), References & File Descriptions
- **Pages 24–27**: Supplementary Tables 1–6 (no figures; all quantitative tables already transcribed in `data.yaml`).
- **Page 29 — SI Notes Fig. S1**: Left out — geometric force-assignment schematic showing binding coordinate $\alpha \in [0, 1]$ and endpoint forces $(1 - \alpha)\mathbf{f}$ and $+\alpha\mathbf{f}$ on filament particles $\mathbf{r}_0^{\mathrm{f}}$ and $\mathbf{r}_1^{\mathrm{f}}$ (no quantitative plot axes).
- **Pages 28, 30–38**: Supplementary Notes text, Table S1 (p. 31), Supplementary References (pp. 34–36), and Description of Additional Supplementary Files (pp. 37–38); no additional figures.

---

## Summary of Key Quantitative Findings from Figure Digitization
1. **Settles `q_gaussian_sf_2pi_convention` & `q_reading_gaussian_sf_printed_vs_table4`**: Evaluating Methods p. 11's literal angular-wavenumber Gaussian (`gaussian_structure_factor_literal_p11`) with Supplementary Table 4's `(A, μ/(2π), σ)` fails (`not reproduced`, up to 100% off-peak error) across all four panels of Supplementary Fig. 7c1–c4, whereas evaluating the same parameters directly in `q/(2π)` spatial-frequency coordinates (`gaussian_structure_factor_freq_variant`, `variant_of: gaussian_structure_factor_literal_p11`) reproduces all four plotted black dotted curves within `0.0015–0.0091` median absolute error (`reproduced` 4/4). In addition, digitizing all 120 blue markers in 7c1–c4 reveals the previously unstated circular binning grid: 30 uniform bin centers `q/(2π) = 0.1, 0.3, ..., 5.9 µm^-1` (`Δ(q/2π) = 0.2 µm^-1`, i.e., `dq = 0.4π µm^-1`).
2. **Settles `q_control_labels_swapped_fig3a`**: Independent digitization of Fig. 3a (`fig3a_controls`), Supplementary Fig. 3b (`figs3b_major_bars`), and Supplementary Fig. 3c (`figs3c_controls`) proves that Supplementary Table 2 and Supplementary Fig. 3b are consistent (`attP40`: `345.6–346.5 nm`, `64.5°`, `161 nm`; `attP2`: `363.4–364.0 nm`, `63.5°`, `166 nm`), whereas both Fig. 3a and Supplementary Fig. 3c swapped the `Control (attP40)` and `Control (attP2)` red labels.
3. **Confirms Table 2's `54.4°` entry in `q_indexing_and_table_crossref_typos`**: Digitizing `fig3a_mutants` places `myosin II^DN` at `(54.42°, 325.6 nm)` vs `RhoA^DN` at `(54.53°, 279.0 nm)`, confirming `54.4°` was plotted directly from the underlying data.
4. **Quantifies `q_fig5f_and_figs5d_significance_claims` item (3)**: Digitizing all 41 time bins of Fig. 4f (`fig4f_coherency`) confirms that `4-M` coherency peaks at `0.2238` (`t = 31 s`) and ends at `0.1891` (`t = 40 s`), never reaching Supplementary Table 3's `0.30` within `t ∈ [0, 40] s`.

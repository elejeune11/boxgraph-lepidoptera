# Forensic Audit Notes — `bashirzadeh2021_source_data`

Complete audit of the public Figshare Source Data deposit (`10.6084/m9.figshare.16337442.v1`, 4 files, 156,518 bytes) accompanying Bashirzadeh et al., *Communications Biology* 4, 1136 (2021) (`bashirzadeh2021`). All calculations are executed live from `$BOXGRAPH_SOURCE` via seven `software: true` Python reader models under `models/` without copying any third-party files into this record.

## File Inventory

| File | Size | Contents | Paper Panels Covered |
| --- | --- | --- | --- |
| `Fig4e_persistence_length_data.xlsx` | 22,092 B | 2 sheets (`Sheet5`, 23 formulas; `P_PL`, 0 formulas) — cell-for-cell identical to `Source_data.xlsx` `Sheet5` & `P_PL` | Fig. 4e, Supp. Fig. 9c |
| `Figs_1_and_2_Pvalues.xlsx` | 27,438 B | 19 sheets (`Sheet1` replicate percentages; `Sheet3`–`Sheet20` Excel Analysis ToolPak Welch t-tests; 0 formulas) | Figs. 1c, 1d, 1f, 2d, 2e, 2f |
| `Source_data.xlsx` | 78,335 B | 7 sheets (`Sheet1`, `Sheet2`, `Sheet1 (2)`, `Sheet3`, `Sheet4`, `P_PL`, `Sheet5`; 1,305 formulas) | Figs. 1b–d, 1f–g, 2b, 2d–f, 3c–d, 4e, Supp. Fig. 9c (+ drafts & unrelated myosin/shape-change scratch tables) |
| `Fig4c-d.mat` | 28,653 B | 15 MATLAB v5 `(N, 2)` float64 arrays (`Values{Actin,Alpha,Fascin}{,2,oval,2oval,Ring}`) | Fig. 4c, Fig. 4d, Supp. Fig. 8e |

## Key Findings

### 1. Fig. 4e & Supp. Fig. 9c Persistence Length (`Fig4e_persistence_length_data.xlsx` & `Source_data.xlsx` `Sheet5`/`P_PL`)
- **Why all three printed Fig. 4e Welch p-values failed in `bashirzadeh2021`**:
  1. **Three buggy Excel SEM formulas in `Sheet5` row 32** inflated the caption's `±` values relative to true `SD/sqrt(N)` (while Fig. 4e's plotted error bars actually show the *true* SEMs):
     - `[0:1:10]` (`A2:A23`, `N=22`, `mean=33.37`, `SD=36.04`): `A32 = A33/SQRT(13) = 9.995` (caption `±10`) vs true `SD/sqrt(22) = 7.683` (+30.1%).
     - `[1:1:10]` (`K2:K19`, **`N=18`**, `mean=54.25`, `SD=44.21`): `K32 = K33/SQRT(15) = 11.414` (caption `±11.4`) vs true `SD/sqrt(18) = 10.420` (+9.5%).
     - `[1:3:10]` (`X2:X27`, `N=26`, `mean=99.71`, full `SD=71.24`): `X33 = STDEV(X2:X25) = 73.610` (omits `X26:X27`) and `X32 = X33/SQRT(24) = 15.026` (caption `±15`) vs true `SD/sqrt(26) = 13.970` (+7.6%).
  2. **Wrong caption sample counts**: `[1:1:10]` has `N = 18` in both `Sheet5` (`K2:K19`) and `P_PL` (`Observations = 18`), not `17`; Supp. Fig. 9c `[3:3:10]` (`AG3:AG17`) has `N = 15` (`95.62 ± 21.86 μm`), not `16`.
  3. **Exact reproduction of `P_PL` and Fig. 4e p-values**: Excel ToolPak's Welch t-test (integer-rounded Welch–Satterthwaite `df`) on raw `Sheet5` reproduces all three `P_PL` outputs to 15 digits: `0.116226` (`df=33`, printed `0.1162`), `0.012547` (`df=42`, printed `0.0125`), and `0.00086695` (`df=38`, **truncated** to `0.0008`). Caption means are likewise truncated (`33.366→33.3`, `39.376→39.3`, `54.251→54.2`).
- **Duplicate bundles & Methods `8 < L < 20 μm` violation**:
  - `[1:3:10]` (`X2:Y27`) duplicates 4 full `(L_p, L)` rows in reverse order (`R4=R22`, `R5=R21`, `R6=R19`, `R7=R18`) plus `X15=X16`, leaving 21 unique bundles (`103.99 ± 15.87 μm`); `[1:1:10]` duplicates `(K3, L3) == (K10, L10)`, leaving 17 unique bundles (`55.61 ± 10.96 μm`). Deduplicated Welch p-values are `0.1305`, `0.0171`, `0.00135`.
  - 25/73 bundles with recorded lengths `L` (34.2%; 20 with `L < 8 μm` down to `4.43 μm`, 5 with `L > 20 μm` up to `52.74 μm`) violate Methods p. 8's `"8 < L < 20 μm"` criterion.
  - Because `L_p` is strongly right-skewed, Mann–Whitney U (`p = 0.01496`) and log-transformed Welch (`p = 0.01619`) are significant for `[0:1:10]` vs `[1:1:10]`, whereas untransformed Welch (`p = 0.1162`) is not.

### 2. Complete Swap of All 18 P-Values Between Fig. 1 and Fig. 2 (`Figs_1_and_2_Pvalues.xlsx`)
- `Figs_1_and_2_Pvalues.xlsx` `Sheet1` matches `Source_data.xlsx` `Sheet2` (Fig. 1) and `Sheet3` (Fig. 2) to `< 5e-6` across all 144 replicate cells, and `Sheet3`–`Sheet20` match Excel's integer-rounded-df Welch t-tests to machine precision.
- **Every single p-value printed on Fig. 1 (`1c, 1d, 1f`, 9 values) actually came from Fig. 2's α-actinin+fascin data (`Sheet3`–`Sheet11`), and every single p-value printed on Fig. 2 (`2d, 2e, 2f`, 9 values) actually came from Fig. 1's α-actinin-only data (`Sheet12`–`Sheet20`)!**
- In addition, `Sheet12` (`0.006397`, true Fig. 1c `0.1`) was misprinted by a factor of 10 as `p = 0.0006` on Fig. 2d, and 4-decimal values were truncated rather than rounded (`0.000472→0.0004`, `0.015065→0.0150`, `0.026163→0.0261`, `0.021749→0.0217`, `0.024161→0.0241`).

### 3. Figs. 1 & 2 GUV Counts & Formula Bugs (`Source_data.xlsx` `Sheet2`, `Sheet3`)
- The Fig. 1c/d (`[389 42 23], [67 45 38], [188 145 128]`) and Fig. 2d/e (`[101 98 29], [147 58 71], [137 82 24]`) brackets are **per-replicate total GUV counts across all diameters** (plus unused extra blocks `J1:L3=67` and `F1:H3=188` in Fig. 1), NOT diameter-bin counts (`[314,70,70]`, `[48,35,25]`, `[155,121,88]` in Fig. 1; `[127,57,44]`, `[176,59,41]`, `[182,33,28]` in Fig. 2). Fig. 1f/2f added the first two replicate totals (`389+42=431` mis-added as `429`). This resolves both `q_fig1f_sample_count_mismatch` and `q_fig1g_vs_fig1f_aster_count_paradox`.
- Fig. 1g (`Sheet2` `M48:N50`) uses denominator `12` (`5/12, 7/12`) for `0.2`, not `14` as printed in the caption.
- `A6` in both `Sheet2` and `Sheet3` omits cell `J11`, shifting Fig. 2b's plotted medium Ring probability from `23/149 = 0.1544` (`T2`) to `17/143 = 0.1189` (`A9`). Both sheets also add a spurious `+1` denominator in `Q43` and `O44`.

### 4. Fig. 3c & 3d (`Source_data.xlsx` `Sheet4`)
- Fig. 3c caption swaps `87` and `80` (`Sheet4` has `80` for `NoFascin>16um` and `87 = 61 + 26` for `fascin`) and misstates the diameter threshold as `>20 μm` when `Sheet4` explicitly analyzes `>16 μm` GUVs. Hardcoded summary cells `K4=34, L4=16` and `J6:L6=(15,7,4)` diverge from the stale `16_25`/`>25` raw table `B4:G6`.

### 5. `Fig4c-d.mat` & Missing Panels
- Contains 15 arrays covering **Fig. 4c** (`ValuesAlpha`, `ValuesFascin`, n=213; plus `ValuesActin` and `*2`, n=209), **Fig. 4d** (`ValuesAlphaoval`, `ValuesFascinoval`, plus unplotted `ValuesActinoval` and `*2oval`, n=360), and **Supp. Fig. 8e** (`Values*Ring`, n=18, undocumented in filename).
- In Fig. 4d, `α-actinin` and `fascin` have nearly identical relative variance (`CV = 3.17%` vs `2.55%`); crosslinker sorting outside the cluster is shown quantitatively via correlation with the unplotted `ValuesActinoval` (`r = +0.421` fascin vs `+0.170` α-actinin; second scan `+0.430` vs `+0.078`).
- **Absent from all Figshare files**: Fig. 3b, Fig. 4g, Fig. 4h (AFINES Binder cumulant), Supp. Figs. 3c, 7b/d, 8b/c, 9b, and 11.

---

## Cross-Record Audit: `bashirzadeh2021_source_data` + `graham2024`

- **Exact pairing trigger**: `pipeline/pairs.py` returned `graham2024 cites Bashirzadeh 2021 (ref-12)` — a **false-positive author-surname/year regex match** (`\bBashirzadeh\b.*\b2021\b`). `graham2024` (Graham et al., *PNAS* **121**(3), e2309152121, 2024) actually cites **Wubshet, Bashirzadeh & Liu (*Mol. Biol. Cell* 32, 1634–1640, 2021)** (twice, as `ref-12` and duplicate `ref-44`), where Yashar Bashirzadeh is second author, and never cites Bashirzadeh et al. (*Commun. Biol.* 4, 1136, 2021) or this Figshare Source Data archive.
- **ModelAtlas relations (`relations.yaml`)**: `0` across all `7 × 39 = 273` model pairs. All 7 `bashirzadeh2021_source_data` models are `role: processing`, `software: true` Excel/MATLAB archive readers over discrete selector indices rather than parametric physical equations.
- **Cited links (`cited.yaml`)**: `0`. All 45 `from_source` datums in `bashirzadeh2021_source_data` already point inline to `bashirzadeh2021`; `graham2024` cites no values from `bashirzadeh2021` or `bashirzadeh2021_source_data` (`~17 μm` and `K_bend = 0.075 pN·μm²` come from Gittes et al. 1993, `ref-42`/`ref-SI2`).
- **Cross-record validations (`validations.yaml`)**: `0` in either direction. None of the 7 spreadsheet/MATLAB readers predict `graham2024`'s 60 `origin: measured` VASP–Arp2/3–VCA droplet datums; conversely, `graham2024`'s 3D CytoSim model is `runnable: false` and lacks α-actinin and fascin (*"Importantly, our system lacks fascin"*, p. 8), and `actin_persistence_length_from_rigidity` (`18.11 μm` at `300 K`) describes single bare F-actin filaments (absent from `Sheet5`), not crosslinked bundles (`33.37–99.71 μm`).

### Quantitative Cross-Study Comparison (`bashirzadeh2021_source_data` vs `graham2024`)

| Observable / Scale | `bashirzadeh2021_source_data` (GUVs, α-actinin ± fascin, `5 μM` actin) | `graham2024` (VASP condensates + Arp2/3 + VCA, `2–5 μM` actin, no fascin) | Key Quantitative Contrast |
| --- | --- | --- | --- |
| Persistence length `L_p` | Crosslinked bundles (`Sheet5`, `N = 109`): `33.37 ± 7.68 μm` (`[0:1:10]`) → `99.71 ± 13.97 μm` (`[1:3:10]`; dedup `103.99 μm`); `[3:3:10]` `95.62 ± 21.86 μm` (`K_eff = 0.138–0.413 pN·μm²` at `300 K`) | Single bare filament (`K_bend = 0.075 pN·μm²`): `18.11 μm` (`300 K`), `17.86 μm` (CytoSim `k_B T = 0.0042`), `~17 μm` text | Bundle means are `1.84×–5.51×` single-filament `18.11 μm` (`N_eff = 1.8–5.5` linear / `1.4–2.4` coupled), yet **50.0% (11/22)** of `[0:1:10]`, **42.9% (6/14)** of `[0:3:10]`, and **21.1% (23/109)** of all bundles fall `< 18.11 μm` (`5.6%–14.3%` with fascin) |
| Confinement diameter `D` & `L_p / D` | GUV bins `7–12`, `12–16`, `>16 μm` (Fig. 3d `>24 μm`, Fig. 4c `49.30 μm`; circumferences `22.0–75.4 μm`). At `16–24 μm`: bare `L_p/D = 0.71–1.06`, bundle `L_p/D = 1.39–6.23` | CytoSim fixed `R_drop = 1 μm` (`D = 2 μm`, `L_max = 2πR = 6.283 μm`, `L_p/D = 8.5–9.1`); expt avg `D ≈ 3 μm` (Fig. 1C weighted mean `2.50 μm` VASP / `4.12 μm` +Arp2/3+VCA, `99%+ < 7.4 μm`) | GUVs are `3.5×–25×` wider; **67/73 (91.8%)** of `Sheet5` bundle contour lengths (`4.43–52.74 μm`, median `10.28 μm`) exceed `graham2024`'s `L_max = 6.283 μm` while remaining shorter than large-GUV circumferences |
| Ring vs aster transition | Fascin stiffening suppresses medium-GUV rings (`50.9%` `Sheet2` → `15.4%` true / `11.9%` `A9` `Sheet3`) and centralizes `>16 μm` asters (`3.75%` central / `62.50%` pole `NoFascin` → `29.89%` central overall, `57.69%` at `>24 μm` `+fascin`) | Arp2/3 shortens median `L_fil` below `D = 2 μm` (`~5.6 → ~0.75 μm`, Fig. 5E), suppressing rings/discs (`34.2%`/`53.6%` → `17.7%`/`7.6%`, `76.1%` shells at `150 nM`, Fig. 4F) at `2:15 μM`, and driving asters (`0% → 76.7%` at `500 nM`, Fig. 6C) at `4:1 μM` | Opposite physical knobs: passive bundle stiffening vs confinement diameter in GUVs vs Arp2/3 filament-length shortening and actin:VASP stoichiometry in `~3 μm` droplets |
| Aster core-vs-spoke sorting | `Fig4c-d.mat`: Fig. 4c central peak-norm α-actinin/fascin = `1.911 / 1.469 = 1.30×`; Fig. 4d circle `CV = 3.17%` vs `2.55%` (`15.88%` actin); exploratory Graham-style spoke/core `E_fascin / E_α = 1.19 / 0.86 = 1.38×` (`1.28×` core window); phase correlation `r(actin, fascin) = +0.421` vs `+0.170` | Fig. 6F spoke enrichment (circular scan `1 μm` outside core): VASP `0.563 ± 0.175` vs Arp2/3 `0.074 ± 0.070` (**`7.6×`**; text `6.8 ± 2.3%` vs `0.69 ± 0.52%`, **`~10×`**) | VASP/Arp2/3 segregation in `graham2024` has `~6×` stronger fold-enrichment amplitude than α-actinin/fascin sorting in `Fig4c-d.mat`, where arm sorting is detected primarily via actin phase correlation |


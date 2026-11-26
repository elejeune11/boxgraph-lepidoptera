# Independent Audit of `tang1996`

Audited `corpora/actin_bundling/records/tang1996/` independently against all 8 page images (`pages/page-001.png`–`pages/page-008.png`) and extracted text (`text/page-001.txt`–`text/page-008.txt`) of Tang & Janmey (1996), *J. Biol. Chem.* 271(15), 8556–8563, plus `corpora/actin_bundling/papers.yaml` for corpus cross-links.

## Transcription Fixes Made

- **`ref-12` (p. 7 / 8562):** Restored printed text `"Manning, G. S. (1978) Q. Rev. Biophys. II, 179-246"` in `references.yaml`. The page image (`page-007.png`) clearly prints Roman numeral `II` (a journal misprint for volume `11`), which the first reader had silently corrected to `11`.
- **`ref-30` (pp. 7–8 / 8562–8563):** Restored printed diacritic `"Käs"` (was `"Kas"`) in `references.yaml` to match `page-007.png`.
- **`chemical_binding_eq4` (p. 7 / 8562):** Updated the `equation` string in `models/chemical_binding_eq4/model.yaml` to `"[Lys_18]_o = [Lys_18]_f + v * [Lys_18]_b = [Lys_18]_f + k * v * [Lys_18]_f"` so it includes the intermediate equality printed in Appendix Eq. 4 on `page-007.png` (no change to `model.py`).

## Disagreements with First Reader (Logged as `reading` Questions; Original Entries Preserved)

- **`alpha_mn_manning_pred` / `q_alpha_mn_manning_range_inferred` (p. 4 / 8559):** The paper never literally prints `"1 to 2"` for Mn²⁺; it prints the general statement `"a slope varying between 1 and N in the log-log plot"` and explicitly substitutes only `N = 3` (`"between 1 and 3"`) for Co(NH₃)₆³⁺. Preserved `alpha_mn_manning_pred` as entered and added `reading` question `q_alpha_mn_manning_range_inferred`.
- **`cmp_eq4_vs_fig3b_6mgml` / `q_eq4_vs_fig3b_computed_comparison` (pp. 3–4, 7 / 8558–8559, 8562):** `cmp_eq4_vs_fig3b_6mgml` places a computed run output (`c_lys18_fig3b_6mgml_comp`) in the `printed` slot, so BoxGraph evaluates it with `printed_text: None`, `difference: None`, and `verdict: no verdict`. Preserved the comparison and added `reading` question `q_eq4_vs_fig3b_computed_comparison`.

## Additional Open Questions Added

- **`q_mixed_K_lys18_40pct_unstated_model` (p. 4 / 8559):** Added `unstated_model` question for the claim that `"an additional 40% or so is neutralized by the excess K+ to make up to 80%"` alongside Lys₁₈ (~40%, exact 43.2%). The paper omits Manning's two-variable mixed-counterion equations (`ref-12`); in the limiting-law (`c -> 0`) two-variable theory, total neutralization stays pinned at `1 - 1/zeta = 64.3%` whenever K⁺ remains condensed, while reaching ~77–80% requires the unstated finite-salt two-variable equation at 150 mM KCl (above the `c < 0.1 M` limit cited on p. 4).
- **`q_atp_vs_adp_actin_charge_state` (pp. 2, 4, 7 / 8557, 8559, 8562):** Added `reading` question noting the unreconciled nucleotide/Pi state of F-actin: p. 2 counts tightly bound ATP to obtain `14` excess negative charges (`11 e` net, `~4 e/nm`), which is used throughout p. 4 for KCl-polymerized F-actin, while p. 7 states that Pi release after ATP hydrolysis lowers filament surface charge density by `~9% (1/11)` (to `10 e`, `3.7 e/nm`, `b = 2.70 Å`, `zeta = 2.63`).
- **`q_unstated_experimental_temperature` (pp. 1–2, 4 / 8556–8557, 8559):** Added `unstated_input` question noting that Materials and Methods never reports the experimental measurement temperature; p. 4 directly reuses `lambda_B = 7.1 Å` from the p. 1 `20 °C` water example.

## Verification Summary (No Transcription Errors Found)

- **`data.yaml` (all 62 entries, pp. 1–7):** Verified every printed value, digit, sign, exponent, approximation marker, unit (including all μM/μm OCR traps on pp. 3, 4, 6), page number, `printed_page`, `from_source`, and verbatim quote against the page images and `text/page-NNN.txt`. Zero digit, sign, exponent, or unit errors found.
- **Models & equations (all 25 models, Eqs. 1–5 and Figs. 3b & 5b embedded fits, pp. 1–7):** Verified every equation as printed on the page images (Bjerrum length p. 1, Eq. 1 p. 1, Fig. 3b linear fit p. 3, Eq. 2 & Eq. 3 p. 4, Fig. 5b power-law fits p. 5, Eqs. 4–5 p. 7), every `model.py` numerical implementation, and exact `run()` vs `symbolic()` agreement across all 24 runnable models. `chemical_binding_atp_eq5` is properly marked `runnable: false`.
- **Runs & comparisons (42 runs, 40 comparisons):** Verified all printed results claimed from the paper's models are compared. Confirmed three uncompared run outputs are benign: `run_theta_dna_tri` (`theta_dna_tri_comp = 92.06%`, p. 1 threshold `"reaches 90%, which requires a valence of 3 or higher"` is an inequality above `theta_dna_condensation_threshold = 90%`, while `N = 2` gives `88.1%`), `run_eq4_vs_actin_6mgml` (duplicates `run_eq4_6mgml`), and `run_fig5b_cohex_150mM` (`2.77 mM`, no scalar Co(NH₃)₆³⁺ onset concentration is printed in text; plotted curves await the figure-digitization stage).
- **`added.yaml` (2 entries):** Verified `d_actin_added` (`7 nm`, alternatives `6, 8, 10`) and `M_actin_added` (`42000 g/mol`, alternatives `41800, 43000`) both have clear bases and realistic alternatives, and no hidden unstated constants exist in `data.yaml` or `models/`.
- **`relations.yaml` (2 stated relations) & `refinements.yaml` (empty, all models deterministic closed-form):** Verified both relations pass `boxgraph relations` (`rel_eq4_vs_actin_same_as_fig3b` and `rel_fig5b_scales_as_manning_onset`), and no numerical/stochastic refinements are needed.
- **`references.yaml` (51 entries) & corpus cross-links:** Verified all 51 references against pp. 7–8 (`page-007.png`, `page-008.png`) and confirmed against `corpora/actin_bundling/papers.yaml` that `tang1996` is the earliest paper in the corpus (no cited reference corresponds to another corpus record).

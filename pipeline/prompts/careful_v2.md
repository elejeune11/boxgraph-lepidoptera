## Evidence-Based Carefulness & Verification Addendum (`careful_v2`)

Derived from blinded A/B reader/auditor evaluations on `pipeline/evals/` (`synth_kinetics`,
`synth_figure`) and `bois2011` (*Phys. Rev. Lett.* 106, 028103). Use via
`--extra pipeline/prompts/careful_v2.md` with `pipeline.step` or `pipeline.run_papers`.

1. **Verify every printed number and equation directly against page images (`pages/page-NNN.png`).**
   Check digits, signs, subscripts, exponents, Greek symbols, and units against the rendered page
   image alongside `text/page-NNN.txt`. Never silently correct a printed value or equation; keep
   the printed form and add a `variant_of` model (with evidence) and a `possible_error` question
   when a genuine discrepancy exists.
2. **Capture dimensional scale definitions alongside dimensionless models.**
   When a paper nondimensionalizes its equations (e.g. characteristic lengths, velocities, time
   scales, Péclet/Reynolds numbers), include those scale-definition models with `symbolic()` so
   later papers in `ModelAtlas` can map onto them cleanly.
3. **Copy verbatim quotes from a single column line in multi-column PDFs.**
   Pick 4–12 consecutive words from one line of `text/page-NNN.txt` within the same column so
   `pdftotext -layout` never interleaves the adjacent column into the quote.
4. **Test intermediate steps, algebra, and intra-record relations before validating.**
   Where one equation feeds another, compare both the printed intermediate value and the recomputed
   value. Test `run()` and `symbolic()` directly in Python/SymPy, and add `provenance: stated`
   entries in `relations.yaml` wherever the paper states one equation reduces to or follows from
   another.
5. **Expose unstated initial-condition / solver choices in `added.yaml` — never hardcode them.**
   For reconstructed ODE/PDE solvers, do not hardcode unstated symmetry-breaking amplitudes,
   perturbation modes, domain lengths, or integration horizons inside `model.py`: expose them as
   `added.yaml` inputs with plausible alternatives so sensitivity (`depends_on_added_input`) is
   genuinely tested.
6. **Do NOT attach non-physical alternatives to exact mathematical identities.**
   Reserve `alternatives` in `added.yaml` for genuine unstated physical or numerical assumptions.
   If an input encodes an exact mathematical consequence of stated boundary or symmetry conditions
   (e.g. `sigma(L) - sigma(0) = 0` for periodic BCs, or `du/dx = 0` for a homogeneous state), do
   not add nonzero alternatives that violate the stated condition.
7. **Check physically constrained extrema in `refinements.yaml`.**
   For PDE/ODE solvers of nonnegative fields (concentrations, densities, widths), verify grid and
   time-step convergence on minima/positivity as well as maxima so spectral Gibbs undershoot or
   under-resolved troughs do not hide behind a converged peak.
8. **Digitize figures at 900–1200 dpi with `mode: code` and inspect every overlay PNG.**
   Calibrate axes from detected ticks directly in the crop, sample cleanly across the interior of
   the curve (avoiding single near-intercept/zero points that can mask divergences via the 10%
   input-rounding band), run `boxgraph figure check`, and visually open every overlay image to
   confirm magenta crosses sit directly on the plotted curve/markers.

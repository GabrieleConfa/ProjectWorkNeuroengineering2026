# Zelechowski, Valle, Raspopovic 2020 — Computational model to design neural interfaces for lower-limb sensory neuroprostheses

**Citation.** Zelechowski M, Valle G, Raspopovic S. *A computational model to design neural
interfaces for lower-limb sensory neuroprostheses.* J. NeuroEng. Rehabil. 17:24, 2020.
doi:10.1186/s12984-020-00657-7. File: `Papers/s12984-020-00657-7.pdf`.

**Role in the project.** The **methodological template**: a hybrid FEM + fibre-model pipeline
that evaluates how many fascicles an intraneural electrode (TIME) can selectively stimulate. Our
project is the same question on the median nerve, with the activating function in place of the
NEURON fibre model and MATLAB PDE Toolbox in place of COMSOL.

## Pipeline (Methods, p. 3–5)

1. **Geometry**: two histological sections of the human sciatic nerve (proximal: 37 fascicles;
   distal: 31 fascicles, 58.2 mm²) segmented in ImageJ, rebuilt in MATLAB, **extruded** along the
   nerve axis to 3D in COMSOL. Perineurium = 3% of the fascicle diameter (Grinberg 2008).
2. **Materials** (Table in text, p. 4):

   | Domain | σ (S/m) | Source cited |
   |---|---|---|
   | Epineurium | 0.0826, isotropic | [42, 55] |
   | Endoneurium | 0.571 longitudinal, 0.0826 transversal | [42, 55] |
   | Perineurium | 0.00088 | Raspopovic 2017 [46] |
   | Surrounding saline | 2 | [42, 43, 55] |
   | Electrode shaft (polyimide) | 6.67·10⁻¹⁴ | [52] |

3. **Domain and BCs**: outer cylinder 16 mm diameter, 15.4 mm long on each side of the centre.
   AS = boundary current source (scaled linearly afterwards — the problem is linear).
4. **Electrode**: TIME, 12/16/20/24 AS on both faces of the shaft, the two faces **shifted by
   half the AS pitch** (asymmetric layout, wider spatial coverage). Shaft 20 µm thick, 380–670 µm
   wide, 18 mm (proximal) or 10 mm (distal) long; **AS = circle of 60 µm diameter**.
   FINE (extraneural) also modelled, AS 0.5 × 0.5 mm.
5. **Equation**: Laplace, ∇·(σ∇V_e) = 0, tetrahedral mesh very fine near the AS, coarse elsewhere.
6. **Fibres**: MRG myelinated-axon model in NEURON, 21 nodes, internode L = 100·D. 1, 3 or 5 fibre
   populations per fascicle depending on fascicle size (< 400 µm, 400–800 µm, > 800 µm). Fibre
   density reduced 50× (11,953 → 240 fibres/mm²) keeping the diameter distribution.
7. **Stimulation**: monopolar cathodic, **50 µs** pulse, charge swept 0.5–60 nC (60 steps);
   bipolar configurations by superposition of single-AS solutions.
8. **Outputs**: recruitment per fascicle; selectivity indices (spatial Sel_i > 0.6 and functional
   Sel_s > 0.9); threshold = charge recruiting **10% of the fibres** of the fascicle.

## Results worth reusing

- Fibres beyond **≈ 2 mm** from a TIME AS were never recruited in the 0.5–60 nC range → a
  natural cut-off radius to save computation (80% time reduction).
- TIME: 20 AS optimal; one TIME reaches ≈ 20–24% of fascicles selectively in a large nerve;
  2 TIMEs ≈ 38%, 3 TIMEs ≈ 54%, 4 no significant gain. Bipolar adds ≈ 9–12%.
- **Validation**: simulated thresholds 6.86 ± 2.07 nC (proximal) and 6.37 ± 2.37 nC (distal) vs
  measured in a trans-radial amputee with TIMEs: **median nerve 5.39 ± 0.98 nC**, ulnar
  6.46 ± 0.72 nC (Petrini et al. 2019). Thresholds were insensitive to nerve geometry because they
  are set at the intrafascicular level. → Direct benchmark for our phase-5 thresholds.
- **[CHECK] Fig. 7 contradicts the text** (seen 2026-10-19): per its legend the light-green bar
  (≈ 5.4 nC) is the *ulnar* nerve and the dark-green one (≈ 6.5 nC) the *median*; the same swap
  holds for proximal/distal sciatic. Verify in Petrini 2019 before quoting 5.39 nC as "median";
  the target range is 5.4–6.5 nC either way (≈ 110–130 µA at 50 µs).
- Other text/figure mismatches: distal TIME-12 is ≈ 17% in Fig. 3b (text: 19.35%, which is the
  TIME-16 bar); distal monopolar FINE-12 ≈ 10% in Fig. 6b (text: 12.58%). The proximal section
  is ≈ 18.7 × 8.8 mm (Fig. 2), not "187 × 88 mm" as written.

## How we use it

- Structure of phases 1–5 (2D section → extrusion → FEM → fibre-level criterion).
- Starting conductivity set (compare with Horn 2023 in `docs/parameters.md`).
- Outer-domain size, AS size, pulse width, 10% threshold criterion, 2 mm range.
- Selectivity indices if we want to go beyond "is every fascicle reached".

## Caveats / differences from our project

- Sciatic nerve (huge, 37 fascicles) vs median (≈ 12–22 smaller fascicles): less tissue to cover,
  one electrode may suffice — that is exactly what we must test.
- Their criterion is a full fibre model (NEURON/MRG); ours is the activating function, a cruder
  but cheaper criterion — note this in the limitations.
- Injected current used for the FEM solves: **2 µA for TIME, 20 µA for FINE** (confirmed with a
  UTF-8 re-extraction of p. 4; the default extraction drops the µ). Irrelevant thanks to linearity.
- Perineurium conductivity is taken from a previous model paper, not from a measurement on human
  tissue; the endoneurium values trace back to older (amphibian/animal) data.

## Newer literature (checked 2026-10-21)

No erratum found. Summarised in `Papers/main.tex` (Zelechowski chapter, "What has been published
since 2020"): `romeni_2023.md` (TIME optimisation in the human median nerve; grounded saline
cylinder r = 10 mm → option for Q5), `strauss_2019.md` (human TIME perceptual thresholds
≈ 37.7 nC per site, ≈ 7 × the 5.39 nC target → issue N3). Also cited there: ASCENT (Musselman
2021), Musselman 2023 (validated VNS thresholds), Marshall 2026 (extrusion vs true-3D: large-fibre
thresholds within median 4.7%, max 19.7%, cuff on vagus), PyFibers (Marshall 2025, 11 fibre models
incl. MRG in NEURON-Python: could check our activating-function threshold), Peterson 2011
(activating-function-based predictor: candidate for Q6), surrogates (Hussain 2024, Toni 2025),
Valle 2022 (thresholds change over months; tissue reaction raises thresholds in models),
Ciotti 2023, Xie 2023, Verardo 2025.

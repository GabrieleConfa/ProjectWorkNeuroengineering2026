# Model parameters and sources

Single reference table for every number used in the model. The values actually used by the code
are set in `config/`; this file says **where each comes from and why it was chosen**. Status:
`candidate` (from literature, not yet chosen) → `chosen` (logged in `decisions.md`).

## Electrical properties (phase 2)

| Domain | Davis 2023 Table 1 (used in the tutor's script) | Zelechowski 2020 (model) | Horn 2023 (measured, canine vagus, 22 °C) | Chosen | Status |
|---|---|---|---|---|---|
| Endoneurium σ_long | 0.571 S/m (ρ = 1.75 Ω·m) | 0.571 S/m | ≈ 1.03 S/m (ρ = 0.971 Ω·m, whole-nerve diameter) | | candidate |
| Endoneurium σ_trans | 0.167 S/m (ρ = 6 Ω·m) | 0.0826 S/m | — | | candidate |
| Perineurium σ | 8.7·10⁻⁴ S/m (ρ = 1149 Ω·m) | 0.00088 S/m | ≈ 2.7·10⁻⁴ S/m (ρ = 3751 Ω·m) | | candidate |
| Perineurium thickness t | tutor: 40 µm for all fascicles | 3% of fascicle diameter (Grinberg 2008) | 13.8 µm on 1.3 mm fascicles (≈ 1%) | | candidate |
| Perineurium sheet resistance ρ·t | ≈ 460 Ω·cm² with the tutor's 40 µm [computed here: 1149 Ω·m × 40·10⁻⁶ m = 0.046 Ω·m²] — close to Horn's measurement | 341 Ω·cm² for a 1 mm fascicle; 133 Ω·cm² for 0.39 mm (median "A" fascicle) | 510 ± 128 Ω·cm² | | candidate |
| Epineurium σ | 0.159 S/m (ρ = 6.3 Ω·m) | 0.0826 S/m | ≈ 0.018 S/m (ρ = 55 Ω·m) | | candidate |
| Surrounding medium σ | muscle 0.086 (trans) / 0.35 (long) S/m (Gielen 1984 via Davis) | 2 S/m (saline) | — | | candidate — saline vs. surrounding tissue to decide |
| Electrode shaft (polyimide) σ | — | 6.67·10⁻¹⁴ S/m | — | | candidate (or model as insulating boundary) |

The values in `config/default_params.m` are, for now, the tutor's (Davis 2023): decision D-c in
`docs/meetings/2026-11-10_kickoff.md`.

Notes
- ρ·t for Zelechowski computed here as (1/σ)·(0.03·D); it depends on fascicle size, Horn's is a
  measured constant. Worth showing in the report.
- Temperature: Horn measured at 22 °C; body is 37 °C. Conductivities increase with temperature.

## Geometry (phase 1) — human median nerve, Delgado-Martínez 2016

| Quantity | Upper arm | Proximal forearm | Distal forearm | Source |
|---|---|---|---|---|
| Nerve cross-section ("epineurium area") | 10.40 ± 1.75 mm² | 7.13 ± 0.59 mm² | 5.80 ± 0.10 mm² | p. 6 |
| Section eccentricity | 0.65 ± 0.02 (all levels) | | | p. 6 |
| Total fascicle area | 4.13 ± 0.23 mm² | 2.46–2.73 mm² | | p. 6 |
| Number of fascicles | ≈ 12 (12.81 ± 0.73, 11.81 ± 0.32) | ≈ 22 (21.87 ± 0.58, forearm) | | p. 7 |
| Fascicle areas | A ≈ 0.13, B ≈ 0.52, C ≈ 1.5–2.0 mm² | almost all A | | Table 1 |

All ± are SEM over blocks of a **single** nerve.

## Upper-arm generator and 3D model (set in `config/pw6_3d_params.m`)

| Quantity | Value | Source | Status |
|---|---|---|---|
| Nerve area (equivalent circle R = 1.711 mm) | 9.20 mm² | Pusnik 2025 (PW6 group document) | chosen (D-09) |
| Nerve shape: ellipse, roundness (minor/major) | 0.62 → semi-axes 2.173 × 1.347 mm | Pusnik 2025 Table 1 (histology, mid-arm), as cited on Marco's branch — **to verify, paper not in `Papers/`** | chosen (D-10) |
| Orientation of the nerve oval in the electrode frame | major axis along x (0°) | assumed — ask tutor | open (D-10) |
| Number of fascicles | 10 | Pusnik 2025 | chosen (D-09) |
| Total fascicle (endoneurium) area | 5.01 mm² | Pusnik 2025 | chosen (D-09) |
| Largest fascicle | 1.49 mm² | Pusnik 2025 / PW6 (Delgado-Martínez cluster C) | chosen (D-09) |
| CV of fascicle diameter | 0.5 | Brill & Tyler 2017: 0.6 ± 0.3 mm, upper arm | chosen (D-09) |
| Inner fascicles (not on the hull) | ≈ 0.58 (N − 6.08) = 2 ± 1 | Brill & Tyler 2017 | chosen (D-09) |
| Gap between fascicles / to nerve surface | 25 µm / 40 µm | assumed (mesh) | to test |
| Perineurium thickness → sheet resistance | 14 µm → ρ·t = 159 Ω·cm² (thin layer) | PW6 / Horn 2023; Zelechowski σ | chosen (D-06) |
| Nerve / cuff length | 70 mm / 50 mm, cuff centred on the contacts | group | chosen (D-02) |
| Cuff wall | insulating void, 0.5 mm thick, shaped on the oval | thickness assumed — ask tutor | open (D-05) |
| Cuff-pad indentation at the face edges | 7 µm → 8.4–13 µm at the face centre (+ local sagitta) | assumed (was 10 µm at the centre on the circle) | to test (D-11) |
| Muscle (outer, grounded) radius | 3 × major semi-axis = 6.5 mm | as 2D (3R), adapted to the oval | chosen (D-10) |
| Contact conductivity (numerical) / idle contacts | 1e3 S/m / 1e-7 S/m (insulating) | numerical stand-ins | to test (D-07) |
| End planes z = ±35 mm | insulated (grounded until 2026-11-27) | D-17: grounded planes caused a fibre-end artefact under an anode | chosen (D-17) |
| Mesh sizes (contacts / fascicles / nerve / outside fine zone / muscle) | 30 / 100 / 150 / 500 / 1200 µm | first guess | convergence study (D-08) |

## Electrode (phase 3)

Provided layout: `data/electrode/sites.mat`, described in `data/electrode/README.md`.

| Quantity | Value | Source | Status |
|---|---|---|---|
| Number of AS | 48 = 8 sections x (4 intraneural + 2 cuff) | tutor, `load_active_sites.m` | given |
| Section positions | z = ±4.3, ±6.3, ±8.3, ±10.3 mm | `sites.mat` | given |
| Intraneural AS radial positions | r = 0.3 and 0.9 mm, both sides of the axis | `sites.mat` | given |
| Cuff AS radial position | R + thickness/2 (face on the nerve surface) | `load_active_sites.m` | given |
| Shank orientations | 0°, 45°, 90°, 135° (sections ∓4, ∓3, ∓2, ∓1) | `sites.mat` | given |
| AS face | 0.2 mm x 0.2 mm | `sites.mat`, slides | given |
| AS thickness | 0.005 mm (5 µm) in `sites.mat`; 50 µm in the slides | `sites.mat` vs PresentazionePW5_1105.pdf | [CHECK] ask tutor |
| AS material | gold | PresentazionePW5_1105.pdf | given |
| Shank body | — | not provided | open |

For comparison, the TIME of Zelechowski 2020: AS circle 60 µm diameter; shaft 20 µm thick,
380–670 µm wide; AS on both faces, staggered by half a pitch.

## Stimulation and activation (phases 4–5)

| Quantity | Value | Source | Status |
|---|---|---|---|
| Pulse width | 50 µs | Zelechowski 2020 | candidate |
| Charge range explored | 0.5–60 nC | Zelechowski 2020 | candidate |
| Internodal length | 100·D | MRG, as in Zelechowski 2020 | candidate |
| Threshold criterion | 10% of fascicle fibres recruited | Zelechowski 2020 | candidate |
| Validation threshold, human median nerve (TIME) | 5.39 ± 0.98 nC | Petrini 2019 via Zelechowski 2020 | reference |
| Activating-function threshold | — | to find | open |
| Fibre diameter distribution | — | to find (median nerve) | open |

## Safety (phase 4–5)

| Quantity | Value | Source | Status |
|---|---|---|---|
| Electrode CIC, Pt | 0.02–0.3 mC/cm² | course slides 08 | reference |
| Electrode CIC, Pt-Ir / TiN | 0.05–0.3 / 0.07–3 mC/cm² [CHECK column alignment in the PDF] | course slides 08 | reference |
| Tissue damage limit | Shannon criterion log D = k − log Q (D µC/cm²/phase, Q µC/phase), k = 1.75 for Pt | Tovbis et al. 2026 citing Shannon 1992 (`literature/tovbis_2026.md`); fitted on macro-electrodes, indicative only for our 0.04 mm² contacts | candidate (verify in Shannon 1992) |
| Anodal block criterion | residual volley at the muscle < 10 % of the target fascicle; other fascicles < 10 % at both ends | D-16 (Tovbis 2026, Zelechowski 2020) | assumed (ask tutor) |
| AS material of the provided electrode | gold | PresentazionePW5_1105.pdf | given (CIC of gold still to find) |
| Electroporation threshold of nerve tissue (field) | — | to find (`p.safety.electroporation_Vpm`; used by `scripts/run_field_safety.m`) | open |
| Field levels reported by the exposure check | 100, 300, 1000, 3000 V/cm; critical volume 1e-3 mm³ | assumed (`p.safety.fieldLevels_Vpm`, `exposedVolume_m3`) | assumed |

## Permittivities (frequency-dependent model, plan: `design.md` roadmap item 5)

| Quantity | Value | Source | Status |
|---|---|---|---|
| Perineurium ε_r | 2018 ± 319 | Horn 2023 | candidate |
| Perineurium sheet capacitance C_s = ε₀ε_r/t | ≈ 0.13 µF/cm² (computed here, t = 13.8 µm); paper reports 13.5 µF/cm² (issue H1) | Horn 2023 | open (H1) |
| Epineurium ε_r | 9.4·10⁶ ± 8.2·10⁶ (Horn); 3.7·10⁶ (Lazorchak, unpublished canine data) | Horn 2023, Lazorchak 2022 | candidate |
| Low/high-frequency ρ, epineurium | 68 / 1.3 Ω·m | Lazorchak 2022 Table 1 | candidate |
| Low/high-frequency ρ, perineurium | 3460 / 273 Ω·m | Lazorchak 2022 Table 1 | candidate |
| Endoneurium ε_r | — | to find | open |
| Muscle ε_r | — | to find (candidate Gabriel 1996, Phys Med Biol 41:2271; not in `Papers/`) | open |

## Gmsh model (D-14; values from Marco's branch, docs/algorithms_and_plan.md)

| Quantity | Value | Source | Status |
|---|---|---|---|
| PDMS lead body cross-section | 0.10 × 0.10 mm | tip 100 µm (tutor, taken as a diameter); width = tip | assumed |
| PDMS sleeve around the nerve | 0.15 mm thick, 50 mm long (= cuff length) | Marco | assumed |
| PDMS conductivity | 1e-7 S/m (stand-in for an insulator, ~1e-12 S/m) | unfilled PDMS (Wacker) | assumed (tutor: filled or not?) |
| Displacement fade along z | 0.4 mm each side of a station | Marco | assumed |
| Cuff pad extra closure | 3 µm beyond full seating | Marco | assumed |
| Edge element size | 5 µm | Marco's edge-refinement test (R converged) | chosen |

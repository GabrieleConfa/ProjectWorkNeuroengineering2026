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
| Tissue damage limit | Shannon criterion (log D = k − log Q) — k and source to verify | — | open |
| AS material of the provided electrode | gold | PresentazionePW5_1105.pdf | given (CIC of gold still to find) |

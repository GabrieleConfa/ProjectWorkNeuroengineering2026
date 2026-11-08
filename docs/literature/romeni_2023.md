# Romeni et al. 2023 — Optimising TIME geometry and insertion in the human median nerve

**Citation.** Romeni S, Losanno E, Koert E, Pierantoni L, Delgado-Martinez I, Navarro X,
Micera S. *Combining biophysical models and machine learning to optimize implant geometry and
stimulation protocol for intraneural electrodes.* J. Neural Eng. 20(4), 2023.
doi:10.1088/1741-2552/ace219. Free preprint: bioRxiv doi:10.1101/2023.02.24.529836 (read on
2026-10-21; the numbers below are from the **preprint**). Not yet in `Papers/`.

**Role in the project.** Closest published work to our question (can TIME active sites reach all
median-nerve fascicles?). Relevant to phases 1, 3, 4 (boundary condition) and 6.

## What they did

- Geometry: histological sections of Delgado-Martínez 2016, "blocks X, XI, XII", above the elbow
  (the 2016 paper only describes blocks I–XI: issue N2 in `Papers/main.tex`). Section extruded
  20 mm.
- FEM (COMSOL): saline cylinder of radius 10 mm, outer surface **grounded**; one solve per
  active site, superposition (lead-field matrix).
- Fibres: MRG-type compartments (node, MYSA, FLUT, STIN); 25,000 fibres per section; activation
  learned by an MLP surrogate.
- Functional topography **assumed**, not measured: 12 muscle groups, distally-exiting groups
  central, neighbouring muscles close.
- Particle-swarm optimisation of the insertion (and inter-site distance) of 1–4 TIMEs; objectives:
  number of fascicles pierced by ≥ 1 AS (M-NPF) or mean distance fascicle–nearest AS (m-AMD).

## Results (preprint)

- Two TIMEs (M-NPF) pierce > 50% of fascicles in the two most proximal blocks; three TIMEs > 50%
  in all blocks and > 75% in the two most proximal.
- Best inter-site distance ≈ 0.5–1 mm; the clinical design (0.75 mm, Valle 2018) loses little.
- Multipolar protocols activate all muscle groups selectively in silico.

## Caveats

- Extruded single section; functional topography synthetic.
- Some symbols (µ, Greek parameter names) are lost in the pdftotext extraction: check the
  injected test current and the electrode parameters in the PDF before quoting them.

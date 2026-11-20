# Davis et al. 2023 — Perineurium thickness, fascicle size and activation thresholds

**Citation.** Davis CJ, Musselman ED, Grill WM, Pelot NA. *Fibers in smaller fascicles have lower
activation thresholds with cuff electrodes due to thinner perineurium and smaller cross-sectional
area.* J. Neural Eng. 20(2), 2023. doi:10.1088/1741-2552/acc42b. Open access: PMC10410695.

**Role in the project.** Perineurium thickness vs fascicle diameter (issue X1), a documented
parameter set (ASCENT defaults) and modelling choices (shrinkage, thin layer, boundary).

## Numbers

Perineurium thickness t vs fascicle diameter d (µm), their Table 2, citing Pelot et al. 2020:

| Nerve | Relation | t for d = 400 µm (our computation) |
|---|---|---|
| Human vagus | t = 0.03702 d + 10.50 | 25.3 µm (6.3%) |
| Human sciatic | t = 0.0275 d | 11.0 µm |
| Pig vagus | t = 0.02547 d + 3.440 | 13.6 µm |
| Rat vagus | t = 0.01292 d + 1.367 | 6.5 µm |

Conductivities (Table 1, S/m): epineurium 1/6.3; endoneurium transverse 1/6, longitudinal 1/1.75;
perineurium 1/1149; saline 1.76; muscle 0.086 (transverse) / 0.35 (longitudinal).

Modelling choices: histology inflated by 25% for shrinkage; perineurium as thin layer (contact
impedance) except multi-bundle fascicles; nerve in a 10 mm-diameter muscle cylinder, outer
boundaries grounded; extrusion 25 mm; threshold convergence < 0.5% for larger domain/mesh.

## Findings

Thicker perineurium → lower peak and broader V_e inside the fascicle → smaller activating
function → higher thresholds; larger endoneurial area lowers current density. Effect size: fibre
diameter > fascicle diameter > fascicle location (cuff).

## Caveats

Human **vagus** and cuff electrodes; no median-nerve relation exists.

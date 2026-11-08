# Decision log

One entry per modelling choice: what was decided, the alternatives, the source, and the
consequence for the results (it feeds the "limitations" part of the discussion).

Format:

```
## D-NN — <short title>  (YYYY-MM-DD, phase N)
Decision: ...
Alternatives: ...
Why / source: ...
Consequence / limitation: ...
```

## Open decisions (to settle before the corresponding phase)

The decisions needed for sprint 1 are laid out with options and a proposal (codes D-a … D-h) in
`docs/meetings/2026-11-10_kickoff.md` §3; once agreed, they get a D-NN entry below. Facts already
checked (2026-11-08, `docs/architecture.md` §5): `femodel` accepts only isotropic σ, `createpde`
accepts a diagonal tensor in 3D; `extrude`, `addCell` and `addVoid` work for the nerve + AS
geometry. The 3D solve is required (stacked 2D solutions are not equivalent, architecture §2).

- **Implantation level** (phase 1): upper arm (≈ 12 fascicles, some large) vs proximal forearm
  (≈ 22 small fascicles). Delgado-Martínez 2016 recommends either for intraneural electrodes.
- **How the cross-section is generated** (phase 1): parametric random generator from the
  Delgado-Martínez statistics vs. segmentation of a histological image.
- **Conductivity set** (phase 2): Zelechowski 2020 vs Horn 2023 vs a mix; which to use as
  baseline and which in the sensitivity analysis.
- **Perineurium representation** (phases 2–3): explicit thin layer in the mesh vs thickened layer
  with rescaled conductivity (same ρ·t) vs thin-layer boundary condition (if the solver allows it).
- **Outer domain and boundary conditions** (phase 3): size, ground vs insulated outer surface,
  location of the return electrode.
- **Electrode geometry in the mesh** (phase 3): explicit shaft volume vs AS-only boundary patches;
  whether PDE Toolbox can build the geometry or an external mesher is needed.
- **Activation criterion and threshold** (phase 5): activating-function form and the threshold
  value per fibre diameter.
- **Number of virtual patients and what varies between them** (phase 6).

## Decisions taken

## D-01 — Location and loader of the provided electrode layout  (2026-11-08, phase 3)
Decision: the tutor's `sites.mat` lives in `data/electrode/sites.mat` (moved from a zip-extracted
folder `src/MatlabCode/sites.mat/`, content unchanged, MD5 in `data/electrode/README.md`). The
project code reads it only through `src/electrode/loadActiveSites.m`, which converts to SI units
and computes the cuff radius as R + thickness/2. The tutor's `load_active_sites.m` is kept
untouched as reference (it points to a non-existent `active_sites.mat`).
Alternatives: fix the path inside the tutor's script and use it as is.
Why / source: project layout (input data in `data/`, SI units inside the code); one loader means
one place to change if the tutor sends a new layout.
Consequence / limitation: none on the results; `tests/testLoadActiveSites.m` fails if a new data
file breaks the stated layout (48 AS, mirror symmetry, radii, orientations).

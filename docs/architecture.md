# Code architecture

**Status: proposal (2026-11-08), to be confirmed by the group at the kick-off meeting
(see `docs/meetings/2026-11-10_kickoff.md`).** Once approved, change the interfaces only by
agreement: other people's code depends on them.

## 1. Principles

1. **One module = one folder in `src/` = one owner** (see `docs/work_plan.md`). The owner
   writes the code, its tests and its entries in `decisions.md` / `parameters.md`; someone else
   reviews it.
2. **Modules talk only through the data structures of section 4.** Field names and units are
   fixed here. If you need a new field, add it here first and tell the others.
3. **SI units inside the code** (m, S/m, A, s, V, C), unit in every variable/field name
   (`radius_m`, `sigmaL_Spm`, `pulseWidth_s`). mm / µm / µA / nC only in plots and tables.
4. **No numbers in `src/`.** Every number comes from `config/default_params.m` (one struct
   `params`) with its source written next to it, and is listed in `docs/parameters.md`.
5. **Linearity:** the volume-conductor problem ∇·(σ∇V) = 0 is linear in the injected current.
   We solve **once per active site with a unit current (1 A)** and scale: V(I) = I · V₁.
   No module ever re-solves the FEM for a different amplitude.
6. **Reproducibility:** every random draw goes through `rng(seed)`; the seed is stored with the
   output. A virtual patient is regenerated from `data/patients/` alone.
7. **Scripts orchestrate, functions compute.** `scripts/run_*.m` only call `src/` functions and
   save to `results/<phase>/`.

## 2. Why a 3D solve (and not stacked 2D slices)

A 2D cross-section model assumes that nothing changes along z: a 0.2 mm active site becomes an
infinitely long strip, and no current flows along the nerve. Solving slices independently and
stacking them would therefore (i) give the wrong spatial decay around a site (2D line source,
V ∝ ln r, instead of a 3D source, V ∝ 1/r), (ii) give V = 0 in every slice without a site,
because the current that flows *along* the nerve — the main path, since the endoneurium
conducts ≈ 3.4× better along the fibres than across them — is missing, and (iii) lose exactly
the quantity phase 5 needs: the activating function is the second difference of V **along z**.

So: the **geometry** is built from a 2D cross-section (extruded along z — slices are fine for
geometry), but the **electrical problem is solved in 3D**. A 2D solve on one section remains
useful as a quick qualitative picture and as a debugging step.

## 3. Pipeline

```
config/default_params.m ──► params
                              │
data/electrode/sites.mat ─► loadActiveSites ──► sites (table, 48 rows)
                              │
seed ─► makePatient ──► patient (nerve radius, fascicles)              [WP1, src/geometry]
            │
            ▼
   buildCrossSection2D ──► xsec (2D geometry + tissue of each face)    [WP1, src/geometry]
            │
            ▼
   buildNerveGeometry3D ──► geom3 (extruded nerve + outer domain)      [WP2, src/mesh]
            │
   addActiveSites(geom3, sites) ──► geom3 with the AS                   [WP2, src/electrode]
            │
   meshNerveModel ──► nerveModel (PDE model + mesh + cell→tissue map)   [WP2, src/mesh]
            │
   tissueConductivity(params) ──► sigma per tissue (σT, σL)             [WP5, src/materials]
            │
   solveUnitCurrent(nerveModel, sigma, siteId) ──► sol (V per 1 A)      [WP3, src/fem]
            │
   makeFibres(patient, params, seed) ──► fibres (positions, diameters)  [WP4, src/activation]
   activatingFunction(sol, fibres) ──► af                               [WP4, src/activation]
   thresholdPerFascicle(af, params) ──► thr (I and Q per fascicle)      [WP4, src/activation]
            │
   checkSafety(thr, sites, params) ──► safety                          [WP5, src/safety]
            │
   scripts/run_pipeline.m ──► results/<phase>/…  (figures, .mat)       [WP5, scripts]
```

Parallel work is possible from day one: WP4 can develop and test the activating function on the
**analytical** potential of a point source in a homogeneous medium, V = I/(4πσr), without
waiting for the FEM; WP3 can develop the solver on a toy geometry (a cylinder with one box)
while WP2 builds the real one.

## 4. Data structures (the interfaces)

Shapes: `(N,1)` column, `(1,3)` row vector. "table" = MATLAB `table`, one row per item.

### `params` — from `config/default_params.m`

| Field | Content |
|---|---|
| `params.geometry` | cross-section generator settings (level, fascicle statistics, perineurium thickness rule, minimum gaps) |
| `params.domain` | outer medium radius, model length along z, outer boundary condition |
| `params.materials` | σT and σL [S/m] of endoneurium, perineurium, epineurium, surrounding medium, electrode/shank |
| `params.electrode` | path of `sites.mat`, which faces inject current, shank yes/no and its size |
| `params.mesh` | target element size per tissue and at the AS [m], growth rate |
| `params.stim` | pulse width [s], current range explored [A] |
| `params.activation` | fibre diameters [m], internode = k·D, nodes per fibre, AF threshold, recruitment criterion |
| `params.safety` | Shannon k, charge-injection capacity of gold [C/m²] |
| `params.patients` | number of virtual patients, list of seeds |

### `patient` — output of `makePatient` (WP1)

| Field | Unit | Meaning |
|---|---|---|
| `id` | – | string, e.g. `"P00"` |
| `seed` | – | seed used (NaN for the tutor's fixed geometry) |
| `nerveRadius_m` | m | nerve radius R (circular nerve; ellipse = open decision) |
| `fascicles` | table | one row per fascicle: `x_m`, `y_m`, `rEndo_m` (endoneurium radius), `tPeri_m` (perineurium thickness) |
| `level` | – | `"upper arm"` / `"proximal forearm"` |

Saved as `data/patients/<id>.mat` together with `params.geometry` used to make it.
`P00` = the tutor's fixed 7-fascicle section (R = 1.5 mm), used until the generator is ready.

### `sites` — output of `loadActiveSites` (done, `src/electrode/loadActiveSites.m`)

Table, 48 rows: `id`, `type` (intraneural/cuff), `section`, `z_m`, `r_m`, `theta_deg`, `phi_deg`,
`thickness_m`, `sideZ_m`, `sideXY_m`, `centre_m` (1,3), `n` (1,3 normal), `w` (1,3 in-plane side).
Full description: `data/electrode/README.md`.

### `nerveModel` — output of `meshNerveModel` (WP2)

| Field | Meaning |
|---|---|
| `pde` | the `PDEModel` object (from `createpde`, see section 5) with geometry and mesh |
| `cellTissue` | (nCells,1) tissue code of each geometry cell: 1 endo, 2 peri, 3 epi, 4 medium, 5 electrode |
| `cellFascicle` | (nCells,1) fascicle index, 0 if not in a fascicle |
| `siteFaces` | (48,1) cell array: geometry face IDs through which each AS injects current |
| `siteArea_m2` | (48,1) area of those faces (to turn the current into a current density) |
| `patient`, `sites` | copies of the inputs, so a saved model is self-contained |

### `sol` — output of `solveUnitCurrent` (WP3)

| Field | Unit | Meaning |
|---|---|---|
| `siteId` | – | AS that injects the current |
| `I_A` | A | injected current (1 A) |
| `result` | – | PDE results object: `interpolateSolution(sol.result, x, y, z)` gives V [V] at any point |
| `V_site_V` | V | mean potential on the AS (for impedance and the safety check) |
| `meshInfo` | – | number of nodes/elements, solve time (for the convergence study) |

### `fibres` — output of `makeFibres` (WP4)

Table, one row per fibre: `fascicle`, `x_m`, `y_m`, `diameter_m`, `internode_m`,
`nodeZ_m` (vector of the z positions of its nodes of Ranvier).

### `af` / `thr` — outputs of `activatingFunction` / `thresholdPerFascicle` (WP4)

- `af`: for each (site, fibre): `afMax_Vpm2perA` = max over nodes of
  (V_{n−1} − 2V_n + V_{n+1})/Δx² for 1 A, and the node where it occurs.
- `thr`: table (site × fascicle): `Ith_A` (current that recruits the chosen fraction of the
  fascicle), `Qth_C` = `Ith_A` · pulse width, `recruitment` curve.

### `safety` — output of `checkSafety` (WP5)

Table, one row per site at its threshold current: `Q_C`, `chargeDensity_Cpm2`, `shannonK`,
`tissueOK`, `electrodeOK` (vs the CIC of gold).

## 5. Technical notes on MATLAB PDE Toolbox (checked 2026-11-08, R2025a, PDE Toolbox 25.1)

- **Anisotropy needs `createpde`, not `femodel`.** The tutor's 2D script uses
  `femodel(AnalysisType="dcConduction")`, whose `ElectricalConductivity` accepts **only a
  scalar** (error: "Value must be scalar or empty"). The endoneurium has σL ≠ σT, so the 3D model
  must use the general form `createpde()` + `specifyCoefficients(..., c=[σT; σT; σL], ...)`,
  one call per cell: a 3-element `c` is accepted in 3D and solves (tested).
- **`extrude(fegeometry(dl), L)`** turns the tutor's 2D `decsg` cross-section into a 3D geometry
  (tested: 1 fascicle → 2 cells).
- **`addCell(geom, box)`** inserts a 0.2 × 0.005 × 0.2 mm box inside a cell of the extruded
  nerve, and the result meshes with `Hmin = 2 µm` (tested on a toy model). This is the way to put
  the intraneural sites in. **`addVoid`** also works, if a site is modelled as a hole with a
  current boundary condition on its faces instead of a gold volume.
- **Not yet tested:** cuff sites touch the nerve surface (face on r = R); `addCell` may refuse a
  cell that touches an existing boundary. A site that crosses a perineurium (as the r = 0.3 mm
  sites do with the tutor's central fascicle) is also untested.
- **Mesh budget:** a 40 µm perineurium meshed at 20 µm along ≈ 25 mm of nerve gives ≈ 10⁷
  tetrahedra — too many. Options are listed in `decisions.md` (perineurium representation).

## 6. Reuse of the tutor's code

`src/MatlabCode/PW_lesson_5_11_students.m` (2D cross-section, tissues, mesh, plots) is the
starting point of `buildCrossSection2D` (sections 4–6 of that script: `decsg`, face → tissue map)
and of the plotting helpers. Keep the original file untouched as a reference; copy the parts
you need into functions.

## 7. Adding code — checklist

- One function per file, `lowerCamelCase`, header: purpose, inputs with units, outputs.
- Validate inputs with an `arguments` block.
- A test in `tests/` for anything with a known answer (analytical case, symmetry, linearity).
- New numbers → `config/default_params.m` + `docs/parameters.md`; new choices → `docs/decisions.md`.
- Run `setup_paths; runtests('tests')` before pushing.

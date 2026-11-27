# The 3D model — what it does, why, and how every piece of code works

Status: **2026-11-13, first version, run end to end on patient UA01, station −4** (results and what
is still unverified in §10).
Decisions referred to as D-NN are in [`decisions.md`](decisions.md); numbers and sources in
[`parameters.md`](parameters.md).

Contents
1. [What the 3D model is, in one page](#1-what-the-3d-model-is-in-one-page)
2. [How we got here — the decisions and the dead ends](#2-how-we-got-here--the-decisions-and-the-dead-ends)
3. [Files and how to run](#3-files-and-how-to-run)
4. [`config/pw6_3d_params.m` — every parameter](#4-configpw6_3d_paramsm--every-parameter)
5. [Virtual patient: `makePatient` + `generateFascicleLayout`](#5-virtual-patient-makepatient--generatefasciclelayout)
6. [Geometry: `buildNerveGeometry3D`](#6-geometry-buildnervegeometry3d)
7. [The finite element method we wrote: assembly, thin layer, interpolation](#7-the-finite-element-method-we-wrote)
8. [The pipeline script `main_pipeline_3d.m`, section by section](#8-the-pipeline-script-main_pipeline_3dm)
9. [Tests](#9-tests)
10. [What is verified, what is not, open questions](#10-what-is-verified-what-is-not-open-questions)

---

## 1. What the 3D model is, in one page

**Question.** For each of the 48 contacts of the electrode, inject a current and find the electric
potential V(x, y, z) everywhere in the nerve. Phase 5 will take the second derivative of V along
each fibre (activating function) to decide which fibres fire.

**Why 3D.** In 2D a 0.2 mm contact is an infinitely long strip, the potential decays like ln r
instead of 1/r, and nothing varies along the nerve, so the activating function (which *is* the
variation along the nerve) cannot be computed. See `architecture.md` §2.

**The physical system.**

```
             z = -35          -25                     0                     +25          +35  [mm]
muscle    |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|  r up to 3a = 6.5 mm
cuff wall               |#############################################|                         nerve oval + 0.5 mm, insulating
nerve     |==================================================================================|  ellipse a x b = 2.173 x 1.347 mm
contacts                          |  8 stations, z = ±4.3 ... ±10.3  |
```

- **Nerve**: an elliptical cylinder (since 2026-11-22, D-10): semi-axes a = 2.173 mm (along x) and
  b = 1.347 mm, area 9.20 mm² (the old circle had R = 1.711 mm), roundness b/a = 0.62; 70 mm long,
  containing 10 fascicles generated from the literature (§5). Same cross-section at every z.
- **Fascicle** = endoneurium (where the fibres are; conducts ≈ 7× better along z than across)
  wrapped in the **perineurium** (14 µm, very resistive: the main barrier for current entering or
  leaving a fascicle). Between fascicles: **epineurium** (connective tissue).
- **Cuff**: a rigid insulating tube shaped on the oval, 50 mm long, tight on the nerve, centred on
  the contacts. Inside it, current cannot leave the nerve sideways — it must travel along the nerve
  and exit at the cuff ends.
- **Contacts**: 48 gold boxes 50 × 200 × 200 µm: per station, 4 on the shank that crosses the nerve
  (r = 0.3 and 0.9 mm, the tutor's positions; the rigid shank is not deflected) and 2 on the inner
  surface of the cuff. The cuff pads are rigid and flat: tangent to the oval and pressed into it
  until the whole face touches the nerve (§2.6, D-11). Only the boxes are modelled, not the shank
  body.
- **Muscle** around everything, up to 3 × the major semi-axis; its outer surface (and the two end
  planes) is the electrical ground (V = 0) — the return electrode is "far away".

**The equation.** At stimulation frequencies the tissues behave as resistors (quasi-static):
charge is conserved and J = −σ∇V, so

  ∇·(σ ∇V) = −f,

with f the injected current per unit volume (non-zero only inside the active contact), σ a
diagonal tensor diag(σ_T, σ_T, σ_L) (z = fibre direction), V = 0 on the ground, no current through
insulating surfaces, and at each fascicle boundary the perineurium condition (§7.2).
**Linearity**: doubling the current doubles V everywhere, so we solve once with 1 A per contact
and scale (control analogy: we compute the "impulse response" of each input channel once).

**What the code does, in order.**

1. Generate the patient's cross-section (`makePatient`).
2. Build a 3D geometry by extruding that section in layers and turning pieces into contacts,
   tissues or holes (`buildNerveGeometry3D`).
3. Mesh it with tetrahedra (PDE Toolbox `generateMesh`).
4. Build the linear system K·V = F ourselves (`assembleConductionP1`, `thinLayerInterface`),
   impose V = 0 on the ground, solve for each contact.
5. Check, plot, save (`main_pipeline_3d.m`).

---

## 2. How we got here — the decisions and the dead ends

This section is the "why" behind the code. Each item is a D-NN entry of `decisions.md`.

### 2.1 What the group decided (2026-11-13)

| Group answer | How it is implemented |
|---|---|
| Nerve 7 cm, cuff 5 cm | `p.domain.nerveLength_m = 70 mm`, `p.cuff.length_m = 50 mm` (D-02) |
| Contacts centred | device centre at z = 0, cuff from −25 to +25 mm (D-02) |
| Cuff tight on the nerve | cuff inner radius = R, no saline gap (D-05) |
| "The electrode cannot bend, so the nerve bends in the vicinity of it" | interpreted as: the electrode is rigid and the nerve takes its shape. (a) inside the cuff the nerve is the circle of the cuff → circular nerve; (b) the flat cuff contacts are pressed 10 µm into the nerve so the whole face touches it (D-03, D-07). **To confirm** that this is what was meant. |
| Points 3 and 4 (ends, perineurium): "do as you want" | ends grounded, to compare with insulated (D-04); perineurium as a thin layer (D-06, see 2.2) |
| Points 5 and 6: proposal (anisotropy, equipotential contacts) | anisotropic σ per element; active contact = equipotential conductor; idle contacts insulating (D-07) |
| Contacts only for now | no shank body (D-07) |
| Everything as parameters | every choice is a field of `config/pw6_3d_params.m` (§4) |
| Upper arm, 10 fascicles, not symmetric, "a biology component" | literature-based random generator (D-09, §5) |

### 2.2 Extrude, not slices + interpolation (D-03)

Asked explicitly. Extrusion was chosen because (i) we have no data on how a given fascicle changes
along 7 cm, so interpolating between slices would invent it, and fascicles that split or merge
cannot be interpolated at all; (ii) the PDE Toolbox has no "loft" operation; (iii) straight
fascicles give straight fibres, which phase 5 needs; (iv) it reuses the 2D section.
Limitation for the report: real fascicles split/merge every few mm.

A **layered** extrusion (`extrude(g2, [h1 h2 ...])`, checked in R2025a) is what makes a 5 cm cuff
on a 7 cm nerve possible with one section: the cuff ring is drawn in the section, and in the layers
inside the cuff it is deleted (void), in the layers outside it becomes muscle.

### 2.3 Dead end 1: meshing the perineurium (→ thin layer, D-06)

First plan: draw the perineurium as a 40 µm ring (instead of the real 14 µm) and scale its
conductivity to keep the same sheet resistance ρ·t. Test (2026-11-13, patient seed 1, 70 mm,
element size 40 µm on the ring): **27.7 million tetrahedra, 4.6 million nodes**, 7 minutes just to
mesh. A thin layer forces small elements in every direction (tetrahedra are roughly isotropic),
along the whole 70 mm, so coarse zones do not help. Not solvable on a laptop.

Solution: the **thin-layer approximation** — the perineurium has no volume in the mesh; the
fascicle boundary is a surface across which the potential jumps (§7.2). It is the "contact
impedance" used by ASCENT (Davis 2023), and Pelot 2019 found it gives thresholds similar to a meshed
perineurium. The PDE Toolbox cannot do this, so we assemble the FE system ourselves (it is ~40 lines
for linear tetrahedra, and tested against exact solutions).

### 2.3b Dead end 3: contacts drawn in the 2D section (→ 3D Boolean union, D-03)

First geometry builder: the contact footprints were drawn in the 2D section and the extrusion had
a layer plane at both faces of every station (21 layers). Every footprint then existed as a cell
in every layer and had to be merged back into its tissue wherever it was not a contact: ≈ 1000
pairwise `mergeCells` on a geometry of ≈ 1500 cells. After 40 minutes it was still running.
Current builder: extrude the tissues only (5 layers), then insert each contact box in 3D with
`union(g, box, KeepBoundaries=[true true])` — the box is split along the tissue boundaries it
crosses, and only its 1–3 pieces are merged (≈ 100 merges in total). Tested on a small case: a box
across a fascicle boundary gives 2 pieces, a cuff box gives the cap in the nerve + the part
embedded in the (deleted) cuff wall, 0.3–1.3 s per box.

### 2.4 Dead end 2: contacts never cut by a perineurium (→ "pierce", D-09)

First generator rule: since the electrode is rigid, the tissue is pushed aside, so no contact is
ever cut by a perineurium. It failed on every seed. Reason: the 16 distinct shank footprints form
two rings of 8 "spokes" (r = 0.2–0.4 mm and 0.8–1.0 mm, every 45°). The fascicles are the same along
z, so all 16 footprints constrain the **same** section (in reality each one deforms the nerve only
locally, at its own z). An average fascicle (r ≈ 0.35 mm) fits neither between two outer spokes nor
outside them, so the only valid layouts have the big fascicle on the axis and the others centred on
spokes: an electrode-shaped, nearly symmetric layout where every contact lands in a fascicle — the
opposite of what we want to test.

Current rule (`contactRule = "pierce"`): the rigid shank goes straight through the tissue, as a
transversal intraneural electrode is implanted. A contact can be in a fascicle, in the epineurium
or across a fascicle boundary (there it cuts the perineurium and the metal shorts it locally).
Only nearly tangent configurations are redrawn, because they would create sliver elements. The old
rule is kept as `"pushAside"` for comparison. **Question for the tutor** (in `project_brief.md`).

### 2.5 Other choices worth knowing

- **Active contact = conductor with a volume source** (D-07). A cell with σ = 1000 S/m and a source
  I/volume spread uniformly inside it. The source integrates to exactly I (so no scaling is
  needed), and the high σ makes the contact equipotential like real metal. Real gold
  (4.1·10⁷ S/m) is avoided because the matrix would become ill-conditioned (entries spanning
  ~10¹¹).
- **Idle contacts insulating** (σ = 1e-7 S/m), as in 2D; `"floating"` is available as an option.
- **Cuff wall = void**: deleting the cells gives an exact insulator and no elements.
- **Ground**: outer muscle surface; end planes insulated since 2026-11-27 (D-17; grounded before, D-04).

### 2.6 The oval nerve and the rigid cuff pads (2026-11-22, D-10, D-11)

Marco's branch `cuff-oval-viewer-convergence` built a second model with an oval nerve and cuff pads
fitted to it. The group took two ideas from it into this model, without its separate pipeline:

- **Oval section.** The median nerve at mid-arm is not round: histology gives a roundness (minor /
  major axis) of 0.62 (Pusnik 2025 Table 1 as cited by Marco — to verify). Same area as before,
  so the semi-axes are R/√0.62 = 2.173 mm and R·√0.62 = 1.347 mm. The fascicle generator packs the
  fascicles inside the oval with the same gaps (distance from the elliptical boundary, §5), so the
  virtual patients and their seeds are kept. `nerveShape = "circle"` gives back the old model.
- **Rigid flat cuff pads; the nerve conforms.** A flat pad cannot touch a curved nerve with its
  whole face unless the nerve is pressed in. `placeCuffContacts` puts each pad at its angle θ from
  `sites.mat`, with the face tangent to the oval there, and pushes it in along that normal until
  both side edges are 7 µm below the native surface. The face centre is then 8.4–13 µm deep
  (deeper where the oval is more curved, at the ends of the major axis). In the geometry the pad box
  simply replaces that thin layer of tissue, so the nerve surface under the pad *is* the flat
  face. On an oval the tangent face is tilted from the radial direction except on the two axes (up
  to ≈ 24° for the pads at 45°): pressing radially, as on the circle, would put the face at an angle
  to the nerve surface.

```
        rigid pad (face flat)                 native oval
      ┌───────────────────────┐            .-'''''''-.
 ─────┴───────────────────────┴─────     .'           '.     the face is tangent to the oval at
        \_ edge 7 µm deep _/             :  nerve        :    the landmark and pushed in along
           centre 8–13 µm                '.             .'    the normal: the nerve conforms,
                                           '-.......-'        the electrode does not bend
```

What was **not** taken from the branch: the shank "deflect" scenario (the shank bends around the
fascicles; here the shank is rigid and its contacts stay at the tutor's positions), the fixed
fascicle template (one patient only), the Gmsh pipeline with smooth shoulders around the pads
(a few µm of difference, a second mesher for the whole group). Taken as helpers: block-wise
assembly (`assembleConductionP1(..., blockSize)`), the pcg convergence check, the fallback path to
`sites.mat`, the regular-grid locator for the viewer export (`locateRegularGrid`).

---

## 3. Files and how to run

| File | Role | New / changed |
|---|---|---|
| `config/pw6_3d_params.m` | all parameters of the 3D model | new |
| `src/geometry/makePatient.m` | virtual patient (was a skeleton) | implemented |
| `src/geometry/generateFascicleLayout.m` | random literature-based section | new |
| `src/mesh/buildNerveGeometry3D.m` | 3D geometry (was a skeleton, new signature) | implemented |
| `src/fem/assembleConductionP1.m` | stiffness matrix, linear tetrahedra | new |
| `src/fem/thinLayerInterface.m` | perineurium thin layer | new |
| `src/fem/interpolateP1.m` | V at any point | new |
| `src/main/main_pipeline_3d.m` | the runnable pipeline | new |
| `scripts/view_model_3d.m` | interactive 3D view of a solved model (geometry + potential planes) | new |
| `src/utils/boundaryFaces.m` | surface triangles of a set of tetrahedra (used by the viewer) | new |
| `scripts/export_viewer_data.m` | exports a solved model for the browser viewer (2026-11-20) | new |
| `viewer/index.html`, `viewer/vendor/` | interactive 3D viewer in the browser, three.js r147 (2026-11-20) | new |
| `tests/testThinLayerFEM.m` | exact-solution tests of the FE code | new |
| `tests/testMakePatient.m` | tests of the generator | new |
| `data/patients/UA01.mat` | patient UA01 (written by the pipeline) | generated |
| `docs/decisions.md` D-02…D-09, `parameters.md`, `architecture.md` "Update 2026-11-13", `project_brief.md` | documentation | updated |

Run (from MATLAB): open `src/main/main_pipeline_3d.m` and press Run, or type
`setup_paths; main_pipeline_3d` in the project root (non-interactive:
`matlab -batch "setup_paths; main_pipeline_3d"`). Time on a laptop: geometry ≈ 3 min, mesh ≈ 1 min,
then one solve per contact.
Tests: `setup_paths; runtests('tests')`.
Change anything (patient seed, lengths, mesh sizes, which contacts to solve...) in
`config/pw6_3d_params.m`, never in the code.

---

## 4. `config/pw6_3d_params.m` — every parameter

The function starts from `pw6_2d_params()` (which itself starts from `default_params()`), so the
materials, the contact size (50 µm thick), the perineurium thickness (14 µm) etc. are the same as
in the 2D model. Tags in the comments: `[source]` literature / PW6 document, `[group]` group choice,
`[assumed]` placeholder or numerical choice to test.

| Field | Value | Meaning |
|---|---|---|
| `patient.id`, `patient.seed` | "UA01", 1 | which virtual patient; the seed makes it reproducible (NaN = tutor's P00) |
| `geometry.generator.nerveArea_m2` | 9.20 mm² | → R = √(A/π) = 1.711 mm (Pusnik 2025) |
| `.nFascicles` | 10 | Pusnik 2025 |
| `.totalFascicleArea_m2` | 5.01 mm² | sum of the endoneurium areas (Pusnik 2025) |
| `.largestFascicleArea_m2` | 1.49 mm² | the one large fascicle (Pusnik / PW6) |
| `.diameterCV` | 0.5 | spread of the other diameters (Brill & Tyler 2017: 0.6 ± 0.3 mm) |
| `.minEndoRadius_m` | 0.12 mm | no fascicle smaller than this (mesh) |
| `.useInnerFascicleRule`, `.innerFascicleTolerance` | true, 1 | 2 ± 1 inner fascicles (Brill & Tyler) |
| `.contactRule` | "pierce" | §2.4 |
| `.contactClearance_m` | 20 µm | minimum distance between a fascicle boundary and a contact edge when they do not cross cleanly |
| `.minGapFascicles_m`, `.minGapNerve_m` | 25 µm, 40 µm | epineurium kept between perineuria / to the nerve surface (mesh; the 40 µm must exceed the 10 µm cuff indentation) |
| `.maxRelaxIterations`, `.maxRestarts` | 4000, 400 | numerical limits of the generator |
| `geometry.layout` | "generated" | |
| `perineurium.representation` | "thinLayer" | §7.2 |
| `perineurium.realThickness_m` | 14 µm | PW6 / Horn 2023 → ρ·t = t/σ = 159 Ω·cm² |
| `domain.nerveLength_m` | 70 mm | |
| `domain.outerRadiusFactor` | 3 | muscle up to 3R |
| `domain.lateralBC`, `domain.endBC` | "ground", "insulated" | end planes: "insulated" since D-17 (was "ground", D-04) |
| `cuff.length_m` | 50 mm | |
| `cuff.wallThickness_m` | 0.5 mm | **assumed**, ask the tutor |
| `cuff.fit` | "tight" | documentation only (the code always uses inner radius = R) |
| `electrode.cuffIndent_m` | 10 µm | flat cuff contact pressed into the nerve (must exceed the 2.9 µm sagitta) |
| `electrode.contactSigma_Spm` | 1e3 S/m | active contact conductivity (numerical stand-in for gold) |
| `electrode.idleContacts` | "insulating" | or "floating" |
| `electrode.insulatorSigma_Spm` | 1e-7 S/m | idle contact conductivity (numerical insulator) |
| `solve3d.sections` | −4 | which stations to solve ([] = all 48 contacts) |
| `solve3d.saveSingle` | true | potentials saved in single precision (half the file size) |
| `solve3d.solver`, `.pcgTol` | "direct", 1e-9 | backslash, or conjugate gradients |
| `mesh3d.fineHalfLength_m` | 12 mm | fine-mesh zone |z| < 12 mm |
| `mesh3d.hContact_m / hEndo_m / hNerve_m / hCoarse_m / hMax_m` | 30 / 100 / 150 / 500 / 1200 µm | target element sizes (§6.7) |
| `mesh3d.hMin_m`, `hGrad` | 5 µm, 1.4 | smallest element, growth rate between neighbours |
| `plot3d.afSpacing_m` | 1 mm | node spacing in the activating-function preview plot |
| `materials.tissueName3d` | names | labels only |

Conductivities (from `pw6_2d_params`, [S/m], order endoneurium, perineurium, epineurium, muscle):
σ_T = [0.0826, 0.00088, 0.0826, 0.13], σ_L = [0.571, 0.00088, 0.0826, 0.56].

---

## 5. Virtual patient: `makePatient` + `generateFascicleLayout`

### 5.1 `makePatient(id, seed, params)`

- `seed = NaN` → returns the tutor's fixed section P00 (unchanged contract).
- Otherwise:
  1. `R = sqrt(area/pi)`, and the contacts are loaded with that R (`loadActiveSites`).
  2. `shankFootprints` (local function) turns the 32 shank contacts into their **16 distinct
     rectangles** in the x-y plane: contacts of sections +k and −k have the same x, y and
     orientation, so duplicates are removed (`unique` on rounded [x y phi]). Each rectangle's corners
     are `centre ± (t/2)·n ± (s/2)·w` with t = 50 µm thickness along the normal n and s = 200 µm along
     the in-plane side w.
  3. `rng(seed, 'twister')` fixes the random stream, then `generateFascicleLayout` is called.
  4. Returns `patient` with `nerveRadius_m`, `fascicles` (table x, y, rEndo, tPeri) and
     `generatorInfo`.

### 5.2 `generateFascicleLayout(gen, peri, footprints)` — the algorithm

Repeated up to `maxRestarts` times until a layout passes all checks:

**(a) Sizes.**
- Draw 9 relative diameters `d = exp(σ_ln · randn)`: a log-normal distribution with median 1.
  For a log-normal, CV² = exp(σ_ln²) − 1, so `σ_ln = sqrt(log(1 + CV²))` gives the wanted CV.
- Areas ∝ d², scaled so that their sum is 5.01 − 1.49 = 3.52 mm² (`a = d.^2/sum(d.^2)*...`).
  Scaling all of them together preserves the relative spread.
- Radii `rEndo = sqrt(area/π)`, the large one first, the others in decreasing order.
- Redraw if any of the nine is ≥ the large one, or smaller than `minEndoRadius_m`.

**(b) Positions: uniform start.** Each centre is uniform in the disc where the fascicle fits:
radius `rMax = R − gap − rOut` (rOut = rEndo + 14 µm), point `rMax·sqrt(rand)` at angle
`2π·rand` (the √ makes the density uniform over the area, not crowded at the centre).

**(c) Relaxation** (at most `maxRelaxIterations` sweeps; stop when a sweep moves nothing):
- *Pairs*: if two fascicles are closer than rOut_i + rOut_j + gap, push them apart along the line
  joining the centres by 1.02× the overlap. The overlap is shared in inverse proportion to the
  areas (`wI = rOut_j²/(rOut_i² + rOut_j²)`), so the big fascicle moves less.
- *Nerve boundary*: a fascicle that sticks out is pulled radially back inside.
- *Contacts*: for each footprint that violates the rule, the fascicle gets a small displacement
  (all the displacements of one sweep are summed and applied together; applying them one by one
  made the moves undo each other). With `"pierce"`: 2 × clearance away from or towards the
  footprint centre, which breaks the near-tangency. With `"pushAside"`: the move that contains the
  footprint, if the footprint centre is inside the fascicle, otherwise the move that excludes it.
- Every 300 sweeps without convergence: a small random kick to all centres (escapes deadlocks).

**(d) Exact final check** (`layoutOK`): gaps between all pairs, inside the nerve, every
fascicle–footprint pair satisfies the contact rule.

**(e) Biology check**: inner fascicles = those not on the convex hull of the centres
(`convhull`); accepted if |nInner − round(0.58(N − 6.08))| ≤ tolerance, i.e. 1–3 for N = 10.

The contact rules (local functions):
- `pierceOK`: the circle of radius r must differ by at least `clr` from the distance to every
  corner and to the closest point of every edge. If the circle passes within 20 µm of an edge or a
  corner without crossing it cleanly, decsg would produce a sliver face.
- `pushAsideOK`: the rectangle is inside r − clr (farthest corner) or outside rOut + clr (closest
  point; `pointRectDistance` returns 0 if the point is inside).

Failures are counted by cause (sizes / relaxation / final check / inner rule), and the error
message prints them if no layout is found.

Measured on seeds 1–10: 1.5–20 s per patient, except seed 7 (147 restarts, 140 s). The total area is
5.010 mm² every time, with 2–3 inner fascicles.

---

## 6. Geometry: `buildNerveGeometry3D(patient, sites, p)`

Input `sites` already has the model thickness (50 µm) and the cuff contacts moved to
r = R − indent + t/2 (done by the pipeline).

### 6.1 Checks on the contacts

- The contacts lie inside the fine zone, which is inside the cuff, which is inside the nerve
  (`max|z_s| + sideZ/2 < zF < Lc/2 < L/2`).
- **Cuff indentation**: δ > (s/2)²/(2R). A flat face of half-width s/2 touching a circle of radius
  R at its middle is (s/2)²/(2R) = 2.9 µm away from the circle at its edges (the sagitta). Pressed
  δ = 10 µm in, the whole face touches the nerve and the nerve circle crosses the two short sides of
  the box at almost 90° (no sliver elements). The cuff contact's centre is at r = R − δ + t/2 (set
  by the pipeline), so its inner face is the chord at R − δ.
- The cuff box fits inside the wall: t − δ < wall thickness.

### 6.2 The 2D section (tissues only)

`decsg` shapes, in this order (the order defines the columns of the table `bt`):
M (muscle circle 3R), C (cuff outer circle R + wall), N (nerve circle R), E1…E10 (fascicle circles
at rEndo — the thin-layer surfaces). The formula "M+C+N+E1+…+E10" unites them, so `decsg` splits the
plane into every elementary piece ("face"), and `bt(f, k) = 1` when face f is inside shape k.
Each face gets `faceTissue`: 1 endoneurium (inside some E_i, and `faceFasc = i`), 3 epineurium
(inside N), 6 cuff-wall ring (inside C only), 4 muscle. The test order matters: the innermost shape
wins. The code asserts that no face lies in two fascicles.

`facePt` (`interiorPoints`): a point safely inside each face, needed later to find the 3D cell that
came from that face at a given z. It is the centroid of the largest triangle of the face in a fine
2D mesh.

### 6.3 Extrusion in 5 layers, cuff wall

Planes: −L/2, −Lc/2, −zF, +zF, +Lc/2, +L/2 → 5 layers with `zone` −2 (outside the cuff, left),
−1 (cuff, not fine), 0 (fine zone, contains all contacts), +1, +2.
`translate(extrude(g2, H), [0 0 -L/2])`: every 2D face becomes one cell per layer. `extrude` starts
at z = 0, hence the translation. The layers exist only to switch the cuff wall on/off and to give the
fine zone its own cells for the mesh sizes.

- In the two end layers, the ring cell is merged with the muscle cell (`mergeGroup`): outside the
  cuff the ring is just muscle.
- In the three cuff layers, the ring cells are deleted (`deleteCell`): void.

### 6.4 Contacts by Boolean union

For each of the 48 sites:
1. `multicuboid(t, s, sideZ)`: a box with x = thickness (50 µm), y = face side (200 µm),
   z = sideZ (200 µm), centred in x-y, from z = 0 to sideZ.
2. `rotate(box, phi, [0 0 0], [0 0 1])`: rotation about the z axis by the angle of the normal, so
   that the local x axis becomes the normal n = (cos φ, sin φ, 0) and y becomes w.
3. `translate(box, centre − [0 0 sideZ/2])`: to its place.
4. `g = union(g, box, KeepBoundaries=[true true])`: the box is added, keeping both the boundaries of
   the existing geometry and those of the box. A box crossing a fascicle boundary therefore comes
   out as one piece per tissue; a cuff box comes out as a thin cap inside the nerve plus a piece in
   the deleted wall region (the contact is embedded in the cuff wall, and its inner face is in the
   nerve).
5. `mergeGroup(g, boxPoints(site))`: the 45 points of a 3 × 5 × 3 grid inside the box (at 80% of
   its half-sides) are located, and all the cells they fall in — the pieces of the box — are merged
   into one contact cell. A contact that crosses a fascicle boundary therefore replaces the
   perineurium there with metal.

`mergeGroup(g, P)`: `mergeCells` takes exactly two cells that share a face and renumbers the cells
after every merge. So the function repeatedly locates the cells containing the points P
(`findCell`), tries pairs until one merge succeeds (a non-adjacent pair throws an error, which is
caught), and repeats until one cell is left. If no pair can be merged, it gives a warning and stops.

### 6.5 Labelling the final cells

- Tissue cells: `facePt` of every face at the mid-z of every layer (except the deleted wall) is
  located with `findCell`, and the cell gets that face's tissue (the ring outside the cuff → muscle),
  fascicle and zone. None of these points can fall in a contact: the layer mid-planes (z = 0,
  ±18.5, ±30 mm) contain no contact.
- Contact cells: the 45 box points must all fall in **one** cell (`assert`). This also checks the
  orientation of the box: if the rotation were wrong, the points along w would leave the box.
- `setRole` (nested function) gives each cell its role and asserts that no cell is reached by two
  different roles; at the end every cell must have a role.

Outputs: `cellTissue` (1 endo, 3 epi, 4 muscle, 5 contact), `cellFascicle`, `cellSite`,
`cellZone`; `footprints` (contact outlines in x-y, used only by the plots).

### 6.6 Faces for the solver and the mesh

- `facesOf(cells)`: union of `cellFaces` of a list of cells.
- **Interface faces** (thin layer) = faces shared by an endoneurium cell and an epineurium cell.
  Faces between endo and a contact, or at the end planes, are not interface.
- **Ground faces**: the outer cylinder, found with `nearestFace` at 8 angles (22.5° + k·45°, away
  from the seams where decsg splits circles) at the mid-z of each layer; plus, if
  `endBC = "ground"`, the faces at z = ±L/2 (found from `facePt` at those z).
- **Element size per face** (`faceSize_m`, smallest wins): nerve faces outside the fine zone
  500 µm; fine epineurium faces 150 µm; fine endoneurium faces 100 µm; contact faces 30 µm; the
  rest (muscle) up to `hMax`.

---

## 7. The finite element method we wrote

The PDE Toolbox builds the mesh; we build and solve the linear system. Reason: the thin layer
(§7.2) needs the nodes of the fascicle boundaries split in two, which the toolbox cannot do.

### 7.1 `assembleConductionP1` — stiffness matrix on linear tetrahedra

Weak form of ∇·(σ∇V) = −f: for every test function φ_i,
∫ ∇φ_i · σ ∇V dΩ = ∫ f φ_i dΩ (+ boundary terms, zero on insulated surfaces).
With V = Σ_j V_j φ_j and φ_j the "hat" functions of linear tetrahedra this gives **K V = F** with
K_ij = Σ_elements vol_e · ∇φ_i · S_e ∇φ_j, S_e = diag(σx, σy, σz) of the element.
Physically K is a conductance matrix (siemens): (K V)_i is the current leaving node i — the same
structure as the nodal analysis of a resistor network.

On a tetrahedron with vertices p1…p4, the hat functions are the barycentric coordinates λ_i,
whose gradients are constant. With a = p2−p1, b = p3−p1, c = p4−p1 and 6·vol = V6 = a·(b×c):

  ∇λ2 = (b×c)/V6,  ∇λ3 = (c×a)/V6,  ∇λ4 = (a×b)/V6,  ∇λ1 = −(∇λ2+∇λ3+∇λ4)

(check: ∇λ2·(p2−p1) = a·(b×c)/V6 = 1, and ∇λ2 is orthogonal to the face p1p3p4).
The code computes these for all elements at once (vectorised), then the 16 entries per element
`vol .* sum((S∇λ_i) .* ∇λ_j)`, and `sparse(I, J, V)` adds up the contributions of elements that
share nodes. An `assert` rejects zero-volume elements.

### 7.2 `thinLayerInterface` — the perineurium as a surface

**Physics.** A layer of thickness t and conductivity σ_p is, for current crossing it, a resistance
per unit area ρ·t = t/σ_p (for us 14 µm / 0.00088 S/m = 0.0159 Ω·m² = 159 Ω·cm²). When t is small
compared with everything else, the layer can be shrunk to a surface with the same resistance per
area: across it the current density is continuous and the potential **jumps**,

  J_n = G_s (V_in − V_out),  G_s = σ_p / t  [S/m²].

Circuit analogy: a sheet of resistors in series between the endoneurium and the epineurium.

**FE implementation.**
1. Every node on the interface (not shared with a contact) gets a **copy**. Elements inside the
   fascicle (endoneurium) are renumbered to use the copies; all other elements keep the originals.
   Now the inside and the outside are two disconnected meshes that touch geometrically.
2. They are connected through the surface term of the weak form, ∫_Γ G_s (V_in − V_out)(φ_in − φ_out) dS,
   which in matrix form is G_s·[M −M; −M M] on (V_out, V_in), with M the surface "mass" matrix of
   the interface triangles. For a linear triangle of area A: M = A/12·[2 1 1; 1 2 1; 1 1 2].
3. **Finding the interface triangles**: take the 4 faces of every inside tetrahedron, sort the node
   triplets, count them (`unique` + `accumarray`). A face counted once is on the boundary of the
   inside region. Keep those whose three nodes are all interface nodes (this drops faces at the end
   planes and faces against a contact).

Where a contact crosses the fascicle boundary, its nodes are excluded from the duplication: the
metal shorts the perineurium there, which is the physics of a contact that pierced it.

Output: the renumbered elements, the extended node list (copies have the coordinates of their
originals), `copyOf`/`dupOf` maps, the interface matrix K_int, the triangles and their areas.
The total system is `assembleConductionP1(tl.nodes_m, tl.tets, sigma) + tl.K`.

### 7.3 `interpolateP1` — V at any point

`pointLocation` on a `triangulation` of the mesh returns the element containing each point and
its barycentric coordinates b; then V(point) = Σ b_k V_(node k). It is exact for linear fields
(tested). Pass the mesh **after** duplication (`tl.tets`, `tl.nodes_m`): a point inside a fascicle
then lies in an endoneurium element, which uses the inside copies, so it gets the inside value.
Points outside the mesh (or in the deleted cuff wall) return NaN.

---

## 8. The pipeline script `main_pipeline_3d.m`

**Section 1 — parameters, patient, contacts.**
- Builds the patient and saves it in `data/patients/<id>.mat` together with `genParams` (generator
  settings, perineurium settings, sites file): that is all you need to regenerate it.
- Prints the fascicle table (position, radius, area).
- Loads the sites with the patient's R, sets the thickness to 50 µm, moves each cuff contact to
  r = R − indent + t/2 (pressed into the nerve).
- `contactLocation`: for each shank contact, 201 points along the mid-line of its face, each
  classified as inside fascicle i or epineurium → e.g. "F3 100%" or "epineurium 35%, F2 65%".

**Section 2 — geometry**: `buildNerveGeometry3D`, prints the number of cells per tissue, faces,
merges, time.

**Section 3 — mesh.**
- `Hface` is built from `faceSize_m`: one {faces, size} pair per distinct size.
- `generateMesh(..., GeometricOrder="linear")`: linear tetrahedra, which is what our assembler
  handles.
- `elemCell`: the cell of every element (`findElements(..., "Cell", c)` for each cell); then
  `elemTissue` and `elemSite` via the cell tables. The code asserts that every element has a cell.

**Section 4 — FE system and solution.**
- `Gs = σ_peri / t_real` (prints the sheet resistance in Ω·cm² as a check: 1e4/Gs).
- Interface nodes = nodes on the interface faces minus the nodes of any contact element; then
  `thinLayerInterface`. Printed check: the interface area should be close to Σ 2π·rEndo × 70 mm
  (a little less where contacts cut it).
- Conductivity of each element: column `elemTissue` of a 3×5 table
  [endo σT σT σL | (unused) | epi | muscle | idle contact]. `K0` = assembly + thin layer.
- Ground nodes: nodes on the ground faces, **plus their inside copies** (at the end planes the
  fascicle boundaries touch the ground, and both copies must be at V = 0). `free` = all the others.
- Contacts to solve: those of `p.solve3d.sections` (default: station −4, 6 contacts).
- For each contact s:
  - if idle contacts are insulating, its own cells switch to the conductor:
    K = K0 + assembly of its elements with σ = (1e3 − 1e-7) S/m (only the difference is added);
  - load vector F: source density I/Vol_s, which for a linear element gives vol_e·f/4 to each of
    its 4 nodes (`accumarray`); Σ F = I exactly;
  - solve K_ff V_f = F_f (V = 0 on the ground) with `solveSPD` (backslash, or `pcg` with an
    incomplete-Cholesky preconditioner);
  - **checks**: current to ground = −Σ_ground (K V − F)_i (the "reaction" at the grounded nodes)
    must equal I (`conservErr_pct`); V spread on the contact (max − min)/mean (`spread_pct`) must be
    small (equipotential); contact potential Vc = mean V on its nodes; for 1 A,
    **R = Vc/I = access resistance** of the contact (printed in kΩ).

**Section 5 — results.**
- First saves `results/main3d/<id>_3d.mat` (see the last item), **before** any figure: a slow or
  interrupted plot must not leave `export_viewer_data` and `export_block_profiles` with an old file.
- `plotSectionLayout`: the cross-section (fascicles labelled F1…F10, all footprints) and a
  z-schematic (nerve, cuff, muscle, stations, layer planes, fine zone).
- **Sampling points located once (2026-11-26).** V is read from the FE solution in two steps: find the
  tetrahedron containing the point (slow), then interpolate linearly with the barycentric weights
  (P1, exact for our elements). Only the first step depends on the geometry, so it is done once:
  `locateOnLines` for the fascicle axes (shared by all 48 contacts) and `locateOnPlane` for each
  station's cross-section, both with `locateRegularGrid` restricted to the elements that touch the
  line or the plane; `evalP1` then applies the stored weights to each contact's V. Same values as
  `interpolateP1`/`pointLocation` (checked: differences ≈ 1e-15 on a test mesh), but ≈ 1 min for the 8
  figures instead of hours (the first full run located every point again for every contact in the
  4.5 M-element mesh, and `pointLocation` cannot be interrupted with Ctrl+C).
- `plotContactResults`, one figure per solved station, one column per contact (each figure is
  closed after saving):
  (1) V on the plane through the contact, for 1 µA, with outlines; title = resistance, subtitle =
  where the contact is;
  (2) V(z) on the axis of every fascicle, ±12 mm around the contact;
  (3) the second difference of V over Δz = 1 mm (one internode of a 10 µm fibre):
  −(V(z−Δ) − 2V(z) + V(z+Δ))/Δ². The minus sign shows it for a **cathodic** pulse, the usual way to
  excite: positive = depolarisation = where a fibre would fire first. This is a *preview* of
  phase 5, not the threshold criterion.
- The saved `results/main3d/<id>_3d.mat` (struct `out`: params, patient, sites, results table, mesh
  after duplication, cell roles, potentials for 1 A) with `-v7.3` (large-file format).

### 8.1 The 3D viewer `scripts/view_model_3d.m`

Reads `results/main3d/<id>_3d.mat` (no recomputation) and opens two rotatable figures, also saved
as `<id>_3d_view_geometry.png` and `<id>_3d_view_AS<k>.png`. Options at the top of the script:
`zRange_mm` (slab to draw — the whole 70 mm × 3.4 mm nerve is a thin stick, so draw the part around
the contacts of interest), `zScale` (compress z), what to show, `showMeshEdges` (to inspect the
mesh), `potentialSite` (which solved contact), `planesOutside`.

How it works:
- The duplicated thin-layer nodes are mapped back to their originals (`Torig`), otherwise the faces
  between a fascicle and the epineurium would not match and every fascicle boundary would look like
  an outer surface.
- Surface of a region = `boundaryFaces` of all its elements (faces that belong to one element only),
  keeping only faces whose three nodes are inside the slab; the cut is marked by the analytic
  outlines of the fascicles, nerve and cuff at the slab ends.
- The cuff wall is not in the mesh (it is a void), so it is drawn as an analytic transparent tube.
- Figure 2 interpolates V (`interpolateP1`, on the mesh with the copies, so values inside a fascicle
  are the inside ones) on the cross-section through the contact and on the longitudinal plane through
  the axis and the contact. By default only the nerve is coloured: inside the cuff zone the muscle
  is separated by the insulating wall and sits near 0 V, which hides the contrast inside the nerve.

What to look for: in the AS 2 view the potential is high along a vertical band that is fascicle F2 —
the current injected inside a fascicle is kept there by the perineurium and spreads along the fibres
(σ_L ≈ 7 σ_T), which is exactly why an intrafascicular contact is selective.

### 8.2 The interactive browser viewer `viewer/index.html` (2026-11-20)

A real-time 3D view for exploring the model and for presenting it: the nerve is a rotatable object
on a black background, fascicles and nerve drawn as bright semi-transparent shells so that the
contacts are visible inside, and the field of a contact switched on with a key or a click as a
glowing see-through volume. The nerve can be cut (dissected) transversally or lengthwise, and the
cut face shows the tissues and the heatmap. Nothing is computed in the browser: it only displays the
FE solution exported once from MATLAB.

**Two steps.**

1. `setup_paths; export_viewer_data` (MATLAB). Reads `results/main3d/<id>_3d.mat` and writes
   `results/viewer/<id>/`:
   - `geometry.js`: the model in mm — nerve radius and length, fascicles (x, y, r), cuff, all 48
     contacts (centre, normal n, side w, sizes, solved or not, location, access resistance), the
     sampling grid;
   - `AS<k>.js`, one per solved contact: **V** (potential), **|E|** = |∇V| (field),
     **|J|** = |σ∇V| (current density) and **AF** (activating function, see below) on a regular
     grid inside the nerve, for 1 µA.
   Options at the top of the script: grid spacing (50 µm across, 100 µm along the nerve), grid
   length (|z| ≤ 12.5 mm, covers all stations), decades of the log scale.
2. Open `viewer/index.html` (double-click; Chrome, Edge or Firefox; works offline).
   `?patient=<id>` loads another exported patient.

**Files the viewer needs** (all in the repository since 2026-11-20, so it works right after
cloning): `viewer/index.html`, `viewer/vendor/` (three.js) and `results/viewer/<id>/` (`geometry.js`
+ one `AS<k>.js` per solved contact, ≈ 7 MB each). `.gitignore` versions these and the 3D figures
(`results/main3d/*.png`) but not the rest of `results/`: the `.mat` solution (50 MB, ≈ 150 MB with 48
contacts: too large for GitHub) and the location cache are regenerated locally.

**How the data are sampled (and why so).**
- Every grid point inside the nerve is located in the mesh *after* the thin-layer node duplication,
  so a point inside a fascicle takes the endoneurium-side value, as in `interpolateP1`. Locating
  ~0.9 M points with `locateRegularGrid` (each nerve or contact element lists the grid points in
  its bounding box; adapted from Marco's branch) takes ≈ 4 s instead of ≈ 13 min with
  `pointLocation`; the result is cached (`gridLocation_cache.mat`, keyed on the mesh, the grid and
  the creation time of the solution). About 0.35% of the points (next to the nerve surface, where
  the meshed polygon lies inside the true ellipse) are not in the mesh and are left empty.
- V is interpolated linearly (P1). **|E| and |J| are taken from the element that contains the
  point**, where the P1 gradient is constant and exact. A finite difference on the grid would
  instead see the potential jump across the perineurium as a huge fake field on every fascicle
  boundary.
- |J| uses the anisotropic σ of the element (σ_z = σ_L in the endoneurium); inside the active
  contact σ = 1000 S/m (D-07), inside idle contacts ≈ 0 (insulating).
- **Activating function (AF)**, the quantity that decides whether a fibre fires: for a myelinated
  fibre along z with nodes every Δz, f = (V(z−Δz) − 2V(z) + V(z+Δz))/Δz². The solution is for a
  positive (anodic) current; for a cathodic pulse, the usual way to excite, the sign flips, so the
  viewer shows AF = −f, where AF > 0 = depolarisation. Δz = `p.plot3d.afSpacing_m` = 1 mm (internode
  of a 10 µm fibre, 100·D), as in the pipeline preview. Computed from the grid V along lines parallel
  to z (fibres are straight in the extruded model), kept **only in the endoneurium** (where the fibres
  are); hyperpolarised points get the lowest colour, the last 1 mm at each grid end has no value. It is
  a map of the stimulus, not the threshold criterion (phase 5 needs fibre diameters and an AF
  threshold). The export also prints, per contact, the AF reached by **10% of the fibres** of each
  fascicle (each grid line along z = one fibre, its maximum along z; 90th percentile per fascicle;
  10% = recruitment level of Zelechowski 2020). The single peak per fascicle is not used: it comes from
  the few fibres touching the contact edges, where the AF is nearly singular (with the peak, AS 1, half
  in F2, would look stronger on F2 than AS 2, fully inside it). Results, patient UA01 in the oval
  (2026-11-22, D-10/D-11), station −4, V/m² per µA:

  | Contact | Most driven fascicle | Next | Ratio |
  |---|---|---|---|
  | AS 1 (F2 100%) | F2 3305 | F3 146 | 23 |
  | AS 2 (F7 100%) | F7 9175 | F8 159 | 58 |
  | AS 3 (epineurium 51%, F1 49%) | F1 1307 | F7 619 | 2.1 |
  | AS 4 (F1 100%) | F1 2577 | F4 133 | 19 |
  | AS 5 (cuff, top) | F8 1023 | F3 507 | 2.0 |
  | AS 6 (cuff, bottom) | F5 1629 | F1 462 | 3.5 |

  Same picture as with the circle: a contact inside a fascicle is very selective (19–58×), a contact
  on a boundary drives two fascicles, a cuff pad drives the superficial fascicles under it but
  spreads more (2–3.5×). The numbers are not comparable one by one with the circle run below: the
  fascicle layout of seed 1 is different in the oval.

  Circular nerve, 2026-11-13 (superseded, kept for comparison):

  | Contact | Most driven fascicle | Next | Ratio |
  |---|---|---|---|
  | AS 1 (epineurium 41%, F2 59%) | F2 1679 | F3 220 | 8 |
  | AS 2 (F2 100%) | F2 2147 | F3 84 | 26 |
  | AS 3 (epineurium 41%, F1 28%, F2 31%) | F2 1109 | F1 924 | 1.2 |
  | AS 4 (F1 100%) | F1 2565 | F7 170 | 15 |
  | AS 5 (cuff) | F4 2866 | F5 724 | 4 |
  | AS 6 (cuff) | F9 2204 | F3 498 | 4.4 |

  The intrafascicular contacts are the most selective; the cuff contacts drive the superficial
  fascicle under them about as strongly per µA, but spread more. Note that V alone hides this: a
  contact inside a fascicle has a *lower* potential (lower access resistance, 8.7 vs 10.0 kΩ), and the
  whole cuffed nerve sits at 1.5–4 mV, which swamps the differences between fascicles.
- Each quantity is stored as an 8-bit code (0 = no data): V linear between its minimum and
  maximum in the nerve; |E|, |J| and AF logarithmic over 3 decades below their 99.9th percentile
  (the top 0.1% of the voxels, right at the contact edges, saturate). The colour bar converts back
  to physical units and scales with the current typed in the panel: the problem is linear, so a
  current of I µA multiplies every value by I. The files are JavaScript, not binary, because a
  page opened from disk may load scripts but not fetch files.

**What the viewer does.**
- Geometry drawn analytically (the 3D model *is* one section extruded along z): cylinders for the
  nerve and fascicles, a tube for the cuff, boxes for the contacts. By default only the data zone
  (|z| ≤ 12.5 mm) is drawn; "full length" shows the whole 70 mm. The muscle is not drawn.
- Heatmap = volume ray marching through a 3D texture: along each pixel's ray, inside the nerve,
  samples are coloured with the colour map and accumulated with an opacity that grows with the
  value, so low values are transparent and the high-value region glows. "Window lo/hi" choose the
  part of the stored range mapped onto the colour bar; "Glow" scales the opacity.
- Dissection: a cut plane removes one side of the shells and of the volume; on the plane itself an
  opaque face shows the tissues (each fascicle's endoneurium in its colour, perineurium as a bright
  outline, epineurium dark, cuff wall grey) and, when stimulating, the heatmap. Transverse
  (slider z, or the station buttons) or longitudinal (angle and offset; by default through the
  nerve axis and the selected contact). Contacts are kept visible unless "keep contacts" is off.
- Selecting a contact moves the rotation centre to its station. The pulse train only modulates the
  brightness to illustrate I(t); with a quasi-static model the field follows the current instantly.

| Input | Action |
|---|---|
| left drag / right drag / wheel | rotate / pan / zoom |
| hover / click a contact | info (station, position, tissue, access resistance) / select and stimulate |
| Space or H | stimulation (heatmap) on/off |
| ← → | previous / next solved contact |
| 1 2 3 4 | V, \|E\|, \|J\|, AF (AF is the default) |
| C, ↑ ↓, F | cycle cut mode, move the cut, flip the side cut away |
| P, M, L, R | pulse train, colour map, full length, reset view |

URL options open a prepared view (for slides or a live demo), e.g.
`index.html?site=2&heat=1&cut=transverse&z=-10.3&flip=1` or
`index.html?site=5&heat=1&q=E&cut=longitudinal` (also: `q=V|E|J|A`, `angle`, `offset`, `cmap=inferno`,
`full=1`, `patient`).

Limits: only the contacts solved by `main_pipeline_3d` have a field (station −4 at the moment); the
file size grows by ≈ 7 MB per contact; the volume rendering needs WebGL2 (any recent browser). The
viewer is a visual aid: numbers for the report come from the MATLAB results, not from screenshots.

---

## 9. Tests

`tests/testThinLayerFEM.m` — exact solutions that linear elements reproduce to round-off:

| Test | Set-up | Exact answer |
|---|---|---|
| `linearFieldAnisotropic` | block 1×1×2 mm, σ = diag(0.1, 0.2, 0.5), 0 V bottom, 1 V top | V = z/L, I = σ_z·A/L |
| `thinLayerJump` | two blocks (σ1 = 0.3, σ2 = 0.1) with a layer G_s = 63 S/m² between | J = 1/(h1/σ1 + 1/G_s + h2/σ2); V_in = J·h1/σ1, V_out = 1 − J·h2/σ2; interface area = 1 mm² |
| `volumeSourceConservation` | source in an inner cube, outer faces grounded | current to ground = I; V ≥ 0 |
| `interpolationExactForLinear` | linear field | interpolation exact; NaN outside |

`tests/testMakePatient.m`: same seed → same layout; areas exactly as in the literature; no overlaps,
inside the nerve; not mirror-symmetric; different seeds differ.

All 9 pass (2026-11-13).

Not part of the repository: the feasibility checks run before writing the code (layered
extrusion, `findCell`/`mergeCells`/`deleteCell`, mesh size of a meshed perineurium = 27.7 M
tetrahedra). Their results are quoted in D-03/D-06.

---

## 10. What is verified, what is not, open questions

Verified: the FE assembly and the thin layer (exact tests), the generator (tests + 10 seeds
inspected visually), and the toolbox operations used by the geometry builder (small cases).

**Oval run, 2026-11-22** (D-10, D-11; patient UA01 = seed 1, new layout in the oval; station −4;
`results/main3d/UA01_3d*.png/.mat`):

| Step | Result |
|---|---|
| Patient | ellipse 2.173 × 1.347 mm, 10 fascicles, 3 inner (target 2 ± 1), 8 restarts |
| Cuff pads | face centre 8.4–13.0 µm into the nerve, edges 7.0 µm; face normal up to 24.0° from radial |
| Geometry | 108 cells (50 endo, 5 epi, 5 muscle, 48 contacts), 537 faces, 30 merges, 222 s |
| Mesh | 863k nodes, 4.56M linear tetrahedra, 134 s (more than the circle: longer perimeter, larger muscle) |
| Thin layer | 140k nodes duplicated; interface area 1606.7 mm² vs 1616.2 mm² expected (−0.6%) |
| Solve | ≈ 3–3.5 min per contact (direct solver) |

| Contact | Where | Access resistance | V spread on contact | Current-conservation error |
|---|---|---|---|---|
| AS 1 shank | F2 100% | 9.95 kΩ | 0.04% | 3e-10 % |
| AS 2 shank | F7 100% | 11.51 kΩ | 0.04% | 5e-11 % |
| AS 3 shank | epineurium 51%, F1 49% | 10.36 kΩ | 0.05% | 4e-11 % |
| AS 4 shank | F1 100% | 9.07 kΩ | 0.04% | 4e-11 % |
| AS 5 cuff (top, flat side of the oval) | nerve surface | 26.91 kΩ | 0.02% | 2e-10 % |
| AS 6 cuff (bottom, flat side) | nerve surface | 29.64 kΩ | 0.02% | 7e-10 % |

The same checks pass (conservation, equipotential contacts, access resistance ≈ 1/(4πσa)). AS 2 has
a higher resistance than the other contacts inside a fascicle because F7 is small (r ≈ 0.27 mm):
its current must cross the resistive perineurium sooner.

**First full run (circular nerve), 2026-11-13** — superseded by the oval run above, kept for
comparison (solution of commit 6d65036):

| Step | Result |
|---|---|
| Geometry | 108 cells (50 endo = 10 fascicles × 5 zones, 5 epi, 5 muscle, 48 contacts), 553 faces, 34 merges, 151 s; every cell got exactly one role |
| Mesh | 510k nodes, 2.84M linear tetrahedra, 45 s (one "poor shape quality" warning from the mesher, not yet located) |
| Thin layer | 155k nodes duplicated; interface area 1557.5 mm² vs 1566.7 mm² expected (fascicle perimeters × 70 mm; −0.6% = polygonal circles + holes cut by contacts) |
| Solve | ≈ 2.5 min per contact (direct solver) |

| Contact | Where | Access resistance | V spread on contact | Current-conservation error |
|---|---|---|---|---|
| AS 1 shank | epineurium 41%, F2 59% | 10.05 kΩ | 0.05% | 6e-11 % |
| AS 2 shank | F2 100% | 8.72 kΩ | 0.04% | 3e-11 % |
| AS 3 shank | epineurium 41%, F1 28%, F2 31% | 8.97 kΩ | 0.04% | 1e-10 % |
| AS 4 shank | F1 100% | 9.03 kΩ | 0.04% | 9e-11 % |
| AS 5 cuff | nerve surface | 37.68 kΩ | 0.02% | 8e-10 % |
| AS 6 cuff | nerve surface | 28.79 kΩ | 0.02% | 1e-9 % |

Sanity checks that pass:
- **Current conservation** to round-off, and the contacts are equipotential (spread < 0.1%).
- **Access resistance**: a contact of size a ≈ 0.1 mm in tissue of σ ≈ 0.1 S/m has ≈ 1/(4πσa) ≈ 8 kΩ;
  the shank contacts give 8.7–10 kΩ. The cuff contacts are 3–4× higher: one face only, and the
  current must travel inside the insulated cuff.
- **Inside the cuff V is linear in z** on every fascicle (figure, middle row). An insulated nerve is
  a 1D resistor: slope = I_branch/(σA) with axial conductance σA = 5.01 mm² × 0.571 + 4.19 mm² ×
  0.0826 = 3.21·10⁻⁶ S·m. The current splits in inverse proportion to the path resistances: the near cuff
  end is 15 mm away, the far one 35 mm, so ≈ 70/30 if only the cuffed length counted, ≈ 65/35 with
  the bare segments beyond (my estimate). With 65/35: 0.20 and 0.11 mV/mm per µA; the FE solution
  shows ≈ 0.2 and ≈ 0.1 mV/mm.
- The two cuff contacts, facing each other across the nerve, differ by 30% (37.7 vs 28.8 kΩ): the
  layout is not symmetric, so the tissue under them differs.

First qualitative reading of the activating-function preview (bottom row; 1 µA cathodic, on the
fascicle *axes* only): the contact inside F2 (AS 2) gives ≈ 1150 V/m² on F2 and < 100 on the others
(selective); contacts across boundaries spread over 2–3 fascicles; each cuff contact strongly
drives the superficial fascicle under it (F4, F9: ≈ 1350 V/m²). Phase 5 will need fibres over the
whole fascicle area, not only the axis.

**Not yet verified:**
- where the "poor shape quality" elements are (likely the 10 µm cuff indentation or the pieces where
  a contact crosses a fascicle boundary) and whether they matter;
- mesh convergence (D-08) — note that the mesh is not identical between runs on the same geometry
  (446k vs 510k nodes in two runs), so the convergence study must show the results do not depend on it;
- ground vs insulated ends (D-04), sensitivity to the cuff wall (D-05) and to the numerical contact
  conductivities (D-07);
- memory: one solve needs many GB (≈ 665k unknowns, direct factorization); the first run was nearly
  stopped by low memory. To do: build the reduced matrix once, save each contact as soon as it is
  solved, try the iterative solver.

**Open questions (for the group / tutor):**
- Is "the electrode cannot bend, so we bend the nerve" what §2.1 says (circular nerve, contacts
  pressed in)?
- Does the shank pierce the fascicles or push them aside (§2.4)?
- Cuff wall thickness and material.
- Shank body size (to add it as an insulating volume).
- Which contacts to solve by default (all 48 is 48 solves).

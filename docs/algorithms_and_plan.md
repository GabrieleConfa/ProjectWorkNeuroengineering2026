# Algorithms and implementation plan (lead insertion, fascicle shapes, PDMS volumes, mesh)

Status of 24 Nov 2026. Written to be learnt: each algorithm has a purpose, inputs, steps, outputs and known limits.
Tags: **[GIVEN]** told by the tutor, **[LIT]** from a paper (see References, full text not always read), **[ASSUMED]** our choice.

## 0. Staged plan and status

| Stage | Content | State (24 Nov 2026) |
|---|---|---|
| A | Mesh: finer elements on the pad edges | **done**, converged at about 5 um, in production for all pads |
| B | Random irregular fascicle shapes | **done and in production** (mesh, solver, figures, viewer) |
| C | PDMS volumes | **done and in production**: PDMS lead bodies and an oval PDMS sleeve around the nerve |
| D | Lead insertion (slide / displace / pierce) | **done and in production**: one random insertion per station, from the entry wall to the opposite wall (anchored end); the Monte Carlo statistics are a separate study |
| E | Figures: fascicles before and after insertion | **done** (study figures and station figures); the 3D model and the viewer use the displaced fascicles |
| F | Full re-run, convergence, report | see the report |

**Production geometry ("insertion" mode, `p.geometry.mode`).** `config/quicklook3DParams.m` builds one reproducible
scenario (`src/insertion/insertionScenario.m`): random irregular fascicles with the area and centre of the circular
layout, and one random lead insertion per station. Figures, the mesh, the solver and the viewer all read this
scenario, so they show the same geometry. `p.geometry.mode = "legacy"` restores the circular fascicles and the
previously approved lead trajectories.

* **Leads.** The lead goes from the entry wall to the opposite wall, where it is anchored (it is not tensioned). The
  path in between comes from the insertion algorithm (section 4). **Both ends stop against the device wrapped around
  the nerve** (the PDMS sleeve, whose inner face is the nerve wall): the lead
  centreline ends exactly on the wall and the PDMS body is cut flush with it, so the lead neither goes beyond nor stops short. After the random walk the anchored end is retracted if it overshoots and
  advanced along the push direction if it stops short; the entry end is brought back to the wall in the same way;
  trials that cannot end against the wall at both ends (for example because a fascicle is in the way) are drawn
  again. Verified by `tests/testInsertionScenario.m` (both ends within 5 um of the wall). In the random walk the tip edge first touches the wall (tip radius 50 um), then the centreline is prolonged to the wall.
* **Fascicles in 3D.** Each fascicle is an extrusion of its outline. Near a station where the lead displaces it, the
  outline moves by the displacement, which fades linearly to zero over `bump_m` = 0.4 mm on both sides of the station
  plane. Undisplaced fascicles are straight prisms.
* **PDMS lead body.** The buffered lead centreline (width = tip diameter, 0.10 mm) extruded over 0.10 mm along z, inside
  the nerve; it replaces the tissue it crosses. The pads keep the geometry of the supplied device file and conduct from
  all their faces (the tutor said that the one-sided emission need not be modelled for now).
* **PDMS sleeve.** A straight oval prism, 0.15 mm thick, around the nerve over the cuff length (50 mm, assumed),
  replacing the circular insulating void used before. Outside the sleeve is muscle.
* **Priorities where solids overlap:** pad > PDMS lead > endoneurium > epineurium > PDMS sleeve > muscle.
* **Insulator.** PDMS is an insulator, modelled with the numerical stand-in sigma = 1e-7 S/m.
* **Pad edges:** 5 um elements on all pad edges (`p.mesh3d.hEdge_m`).

Assumed values (all adjustable): lead cross-section 0.10 x 0.10 mm, sleeve thickness 0.15 mm, sleeve length 50 mm,
displacement fade length 0.4 mm. Known limits: the insertion forces are uncalibrated (section 4); a displaced fascicle
cannot push its neighbours; the simulation is 2D per station (the lead is not tilted along z); Gmsh reports a few
ill-shaped tetrahedra (below 0.01 %), mostly in the thin PDMS layers.

Rule: anatomy, leads and cuff are changed only with the data listed in section 6, and every change must show up
consistently in figures, mesh and viewer (same geometry everywhere).

## 0a. Patients, functions and motor patches (`src/patient`)

The geometry comes from a **patient** (virtual or specific, see `patient_data.md`). Each fascicle has a *composition*: the
fraction of its fibres in each functional group (`functionCatalog`: PT, FCRPL, FDS, AIN, THENAR, PCB, DIGITAL), because at
mid-arm fascicles are mixed and no published map assigns fascicles to muscles (`function_map_evidence.md`). For virtual
patients the composition is drawn by `assignFunctions` under a stated hypothesis (random, forearm rule, MR neurography);
about 13 percent of the fibres are motor (Delgado-Martinez 2016) and lie in 1 to 3 random **patches** inside the fascicle.
Motor patches do not change the conductivity (effective endoneurial value); they tell where the motor fibres run, which the
recruitment calculation will need together with fibre diameters and thresholds (not available yet). Figures colour each
fascicle by its dominant function; the viewer labels it. The reference patient V001 gives exactly the geometry used so far.

## 0b. Sensitivity to model choices (`scripts/run_pdms_sensitivity.m`)

Access resistance of AS1 (internal) and AS5, AS17 (cuff), and where the 1 uA goes (current reaching the grounded nerve
end planes at z = +/-30 mm and the grounded outer muscle surface). The meshes are those of the production geometry.

| Case | R AS1 (kOhm) | R AS5 | R AS17 | end planes | muscle surface |
|---|---|---|---|---|---|
| PDMS insulator (baseline), sleeve 50 mm | 18.0 | 32.7 | 33.1 | 69 % | 31 % |
| PDMS conductive 1e-3 S/m, sleeve 50 mm | 15.4 | 28.7 | 28.2 | 10 % | 90 % |
| PDMS conductive 1e-1 S/m, sleeve 50 mm | 13.1 | 7.7 | 7.4 | 4 % | 96 % |
| PDMS insulator, sleeve 36 mm | 16.7 | 31.7 | 31.9 | 24 % | 76 % |
| PDMS insulator, sleeve 22 mm | 15.2 | 29.5 | 30.7 | 7 % | 93 % |
| PDMS insulator, sleeve 50 mm, nerve end planes insulated | 18.0 | 32.8 | 33.1 | 10 % | 90 % |

Reading: (1) **Whether the PDMS conducts matters a lot**: unfilled PDMS is an insulator (bulk conductivity of the order
of 1e-12 S/m; conductive PDMS needs carbon or similar fillers), but if the tutor's PDMS were filled and conductive
the cuff-pad resistance would fall by 12 to 76 % in these two cases. The model uses the insulator by default
(`p.electrode.pdmsSigma_Spm`). (2) The **sleeve length** matters moderately (10 % lower R at 22 mm). (3) The
**boundary condition on the nerve end planes does not matter for R** (insulated or grounded gives the same R within
0.5 %), although most of the current reaches them in the baseline: the voltage is determined by the spreading near the
pad. (4) With the long insulating sleeve the muscle around the sleeve carries little current: the current leaves the
nerve at the sleeve ends (|z| > 25 mm) or at the end planes. These are the reasons for the tutor questions on the PDMS
type and on the real length of the sleeve.

## 1. Stage A: refining the pad edges (mesh)

**Purpose.** The current of a flat conductor crowds on its edges and corners. A coarse mesh smooths this crowding,
so the pad behaves as if it were wider and the access resistance R is underestimated.

**Method.** `p.mesh3d.hEdge_m` (0 = off) adds a Gmsh size field around the pad edges
(`src/mesh/mesh_personal.py`): a `Distance` field to the edge curves and a `Threshold` field that grows from
`hEdge` to the coarsest size over 0.15 mm. The pads to refine can be selected (`p.mesh3d.hEdgePads`).

**Lesson learnt (Gmsh).** A `Threshold` field returns `SizeMax` *everywhere* beyond `DistMax`. If `SizeMax` is small,
the whole model is refined and the mesh explodes (4 million nodes). `SizeMax` must be the coarsest size.

**Result (pads AS1 internal, AS5 and AS17 cuff, production mesh otherwise).**

| Edge size | Tetrahedra | R AS1 (kOhm) | R AS5 | R AS17 |
|---|---|---|---|---|
| off | 548,832 | 16.54 | 31.46 | 30.62 |
| 10 um | 556,982 | 16.89 | 32.94 | 32.05 |
| 5 um | 563,596 | 16.94 | 33.02 | 31.90 |
| 2.5 um | 569,263 | 16.93 | 33.01 | 31.95 |

Refining the edges raises R, and below about 5 um it stops changing (changes under 0.5 %). The edge explanation is
therefore **supported**. The production mesh without edge refinement underestimated R by about 2 % (AS1), 5 % (AS5)
and 4 % (AS17). Three pads cost about 3 % more tetrahedra, much cheaper than refining the whole model. Production now
uses `p.mesh3d.hEdge_m = 5 um` on all 48 pads. Files: `scripts/run_edge_refinement_test.m`,
`results/convergence/edge_refinement.*`. The convergence study scales the other sizes and keeps the edge size fixed.

## 2. Stage B: random irregular fascicle shapes

**Purpose.** Fascicles are neither circles nor ellipses [GIVEN]. Generate plausible random shapes.

**Target statistics [LIT].** Fascicle roundness 0.68 +/- 0.08 (Pusnik et al. 2025, tables 1 and 4; roundness =
4 A / (pi * major_axis^2)). Areas keep the aggregate constraints (10 fascicles, total 5.01 mm2, largest 1.49 mm2).
Circularity from area and perimeter is also used in cadaver morphometry (Ref. 9).

**Algorithm (radial Fourier shape).**
1. Start from an ellipse with the target area and an axis ratio sampled so that the roundness is about 0.68.
2. Perturb the radius: r(theta) = r0(theta) * (1 + sum_k a_k cos(k theta + phi_k)), k = 2..4, with random
   a_k (decaying like 1/k) and random phases. The perturbation keeps the curve smooth and star-shaped.
3. Rescale so the area is exact. Compute the roundness. Reject the sample if it is outside 0.68 +/- 0.08.
4. Place the shapes without overlap, keeping the perineurium thickness (14 um) and a clearance from the nerve surface
   (as `midArmLayout` does now with circles).
5. Fixed random seed, so every run is reproducible.

**Output.** One closed polygon (about 128 points) per fascicle. The mesher needs a spline or polygon instead of the
four ellipse arcs now used in `mesh_personal.py`; fascicle outlines in figures and in the viewer must use the same
polygon. **Limits.** The shape statistics come from one roundness value; real fascicles also have concavities.

## 3. Stage C: PDMS volumes

**Given [GIVEN].** The leads are PDMS. A thin oval PDMS layer, slightly thicker than the lead, rolls around the nerve.
Muscle surrounds it. The tutor said the lead volume is not required, but we try it.

**Model.**
- Unfilled PDMS is an insulator (Refs. 4, 5). In the code it is a numerical stand-in, sigma = 1e-7 S/m
  (`p.electrode.insulatorSigma_Spm`), as for idle pads; true values are far lower, which makes no difference here.
- **Lead body:** a thin strip (width and thickness to be given) around each internal pad along the lead path. The pads
  are **not completely wrapped**: the active face is exposed to the tissue, the back and sides are covered by PDMS [GIVEN].
- **Sleeve:** the offset of the nerve contour by the sleeve thickness (oval, not a circle), replacing the circular
  void of 0.5 mm used now.
- Parameters to add to `quicklook3DParams.m`: tip diameter 0.100 mm [GIVEN, assumed to be a diameter], lead width and
  thickness, sleeve thickness, sleeve length.

**Consequence to check.** The sleeve insulates the nerve surface, so current must leave along the nerve axis; sleeve
length and the grounded end planes then control R. This is tested by a sensitivity run.

## 4. Stage D: lead insertion simulation

**What happens [GIVEN].** The clinician pushes the lead with forceps (a hand force). It enters on one side and is only
fixed at the end, **not** tensioned. At a fascicle three outcomes are possible: the lead is deflected, the fascicle is
displaced, or, for a very large fascicle, the lead passes through it.

**Inputs.** Tip diameter 0.100 mm [GIVEN]. PDMS Young's modulus about 2 MPa [LIT], range 1.3 to 3 MPa for Sylgard 184
at 10:1, tensile (Ref. 3); the compressive modulus is much higher, so state the test mode. Perineurium puncture
threshold F_p [LIT, see limits]. Fascicle mobility (epineurium stiffness) [to find]. Entry side and angle [random].

**Algorithm (Monte Carlo, per station).**
1. Sample: entry side, entry angle, lateral offset, small initial imperfection of the lead. Fixed seeds.
2. Move the tip forward in small steps along the commanded direction (displacement control: the hand imposes the
   motion; the lead is soft).
3. When the tip meets a perineurium compute the encounter angle and three resistances:
   - `F_pierce(theta, tip)`: force to puncture the perineurium (grows when the angle is small, i.e. the lead slides;
     decreases with a steeper angle [LIT, Ref. 7]).
   - `F_displace`: force to move the fascicle sideways. It grows with fascicle size and with the lack of free room
     (neighbours and nerve boundary), and is limited by the epineurium.
   - `F_deflect`: lateral force to bend the lead, small for PDMS (beam theory, E I / L^3).
4. The smallest force wins. A very large fascicle has a large `F_displace`, so piercing becomes the easiest path.
5. Update the geometry: bend the path, translate the fascicle (area is conserved and overlaps with neighbours are
   resolved), or record the intrafascicular segment (fascicle number, entry and exit point).
6. Place the pads along the path at their fixed distances from the tip (`sites.mat`) and record the tissue of each pad
   (epineurium or the endoneurium of fascicle F#).
7. Repeat N times (500 to 1000 per station) to obtain outcome frequencies.
8. Run the FEM only for representative cases (one per outcome), because each mesh takes minutes.

**Honesty rule.** Until `F_pierce`, `F_displace` and the entry statistics are calibrated, the frequencies are
indicative and must not be quoted as real probabilities. A plausibility check: pig data show about 16 % of the path
intrafascicular and a higher penetration with a larger angle (Ref. 8; preprint).

**Code layout (new).** `src/insertion/`: `sampleInsertion.m`, `advanceTip.m`, `encounterForces.m`,
`resolveEncounter.m`, `placePads.m`, `runInsertionMonteCarlo.m`; one script `scripts/run_insertion_study.m`.

**Results of the prototype (100 insertions per station, parameters in `config/insertionParams.m`).**
Contacts always happen (the section is full of fascicles). With the assumed hand force (median 20 mN, far above the
2 mN puncture force) most head-on contacts pierce, because packed fascicles leave little room to be displaced; grazing
contacts (about one third) slide around. The sensitivity panel shows that the outcome depends strongly on the ratio
hand force / puncture force: with a median force of 1 mN displacement dominates, with about 2 mN the three outcomes
are mixed, and above 5 mN piercing saturates near 60 %. **This is why the frequencies must not be quoted as
probabilities until the forces are calibrated.** Large fascicles (F1, F5) are crossed most often, as the tutor expects.
Limits of the prototype: 2D (one section per station), a displaced fascicle cannot push its neighbours (it is
blocked instead), friction and stiffness are single numbers, the nerve is not deformed by the lead.

## 5. Stage E: figures before and after insertion
Implemented in `scripts/run_insertion_study.m`: `insertion_before_after.png` (one section before and after), `insertion_all_stations.png`
(one insertion per station) and `insertion_statistics.png` (outcomes, fascicles crossed, sensitivity to the hand force).
The mesh and the viewer must read the same post-insertion polygons once stage B/E are integrated (see section 0).

## 6. Data still needed
- Meaning of "tip 100 um" confirmed as a diameter; lead width, thickness, length; sleeve thickness and length (now assumed, see section 0).
- Perineurium puncture force for a 100 um tip, and its dependence on angle (the 2 mN of Ref. 6 is for a tungsten
  needle in rabbit nerve; the cadaver range in Ref. 7b is much larger and not comparable).
- Fascicle mobility: stiffness of the interfascicular epineurium.
- Forceps force range and the entry-angle distribution.
- Cuff closing force, or pad gap, for the cuff stage.

## 7. References (search results of 24 Nov 2026; check full texts before citing numbers)
1. Pusnik et al. 2025, Sci Rep, Table 1 and 4: nerve and fascicle area, roundness (already in the literature pack).
2. Delgado-Martinez et al. 2016, Front Neurosci 10:286, fascicular topography of the human median nerve.
3. PDMS elastic properties: influence of fabrication protocol and test method, arXiv:2404.16960; Sylgard 184 tensile
   modulus 1.3 to 3 MPa depending on cure (Univ. Hertfordshire record 906746).
4. Wacker, "Selected electrical properties of silicone elastomers" (insulating; resistivity measurement caveats).
5. Unfilled PDMS is non-conductive; conductive composites need fillers (Sci Rep 2019, s41598-019-49939-8).
6. Jensen et al. 2001, "Measurement of intrafascicular insertion force" (tungsten needle, rabbit nerve, > 2 mN, dimpling).
7. (a) Sergi, Jensen, Yoshida, "In vivo interactions between tungsten microneedles and peripheral nerves"
   (arXiv:1308.0525): insertion force vs tip angle and minimum puncture diameter. (b) Univ. Dundee Thiel cadaver
   record: epineural penetration 0.6 to 16.8 N, fascicles "spin away" (not comparable with 6).
8. MiiNS preprint 2025 (pig, in the literature pack): about 16 % intrafascicular, angle dependence.
9. Cadaver histomorphometry of median, ulnar and radial nerves, PMC5712902 (fascicle circularity, area).
10. Akouissi et al. 2022, J Neural Eng 19:046017, mechanical interaction nerve and intrafascicular implants.
11. Raspopovic et al. 2011, IEEE TNSRE, FEM model of the TIME electrode in rat sciatic nerve.
12. Main et al. 2011, J Biomech: Ogden parameters of the human median nerve (convex platens, wrist).
13. Sala-Blanch et al. (intraneural needle placement is more often interfascicular, human sciatic).

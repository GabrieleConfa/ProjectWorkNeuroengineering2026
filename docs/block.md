# Blocking the signal towards the amputated limb — conduction-block study

Status: first round (2026-11-25) and second round (anode on the same contacts, 2026-11-26) on the
PDE Toolbox model; third round written 2026-11-27 (travelling volley + anodal block judged per fascicle,
`notebooks/block_volley.ipynb`, §5), to run on the Gmsh solution. Choices and their reasons:
`decisions.md` D-12, D-16. Roadmap item: `design.md` §5.3.

## 1. The clinical problem

After an amputation with **targeted muscle reinnervation (TMR)** the surgeon reconnects the cut
median nerve to a new muscle (e.g. a head of the biceps). The motor fibres of the median nerve then
drive that muscle, and its EMG is used to control the prosthesis.

The electrode in the upper arm stimulates the median nerve to give the patient a **sensation**
(sensory feedback from the prosthesis). But a fibre excited by an electrode fires an action potential
(AP) **in both directions**:

- towards the brain (proximal): the wanted sensation (sensory fibres);
- towards the hand (distal): it reaches the reinnervated muscle (motor fibres, and the antidromic
  branch of the sensory ones). The muscle twitches, which is uncomfortable and corrupts the EMG that
  controls the prosthesis.

The goal is to let the proximal-going APs through and stop the distal-going ones.

## 2. The two sites and the proposed scheme

The 48 contacts are in two mirror-image halves along the nerve, the two **sites**:

| Proximal site (z < 0, shoulder) | Distal site (z > 0, hand) | angles: shank / cuff |
|---|---|---|
| section −4, z = −10.3 mm (contacts 1–6) | section +4, z = +10.3 mm (43–48) | 0/180° · 90/270° |
| section −3, −8.3 mm (7–12) | section +3, +8.3 mm (37–42) | 45/225° · 135/315° |
| section −2, −6.3 mm (13–18) | section +2, +6.3 mm (31–36) | 90/270° · 0/180° |
| section −1, −4.3 mm (19–24) | section +1, +4.3 mm (25–30) | 135/315° · 45/225° |

Each station is a TIME shank with 4 intraneural contacts plus 2 cuff pads. Which end faces the
shoulder is not known; with the mirror symmetry it does not matter, and z < 0 is called proximal.

The scheme, as described by the supervisor (2026-11-25): **the device uses an electrode
configuration selected to excite a particular fascicle, and on the other side the same electrodes
work as cathodes to cancel the signal.** So the proximal site stimulates a pattern of contacts aimed
at one fascicle (the "target"), and the distal site drives **the same pattern** (the mirror contacts,
section −k → +k), again cathodic, so that the distal-going AP is cancelled. Two mechanisms can do
this with a cathode, depending on its amplitude (both simulated, §6):

- at about the stimulation current: **collision** (below);
- at a much higher current: **cathodal block** (below).

**Collision.** Two APs travelling towards each other on the same axon annihilate: behind each AP the
membrane is refractory (Na channels inactivated), so neither can pass the other. This is why
repeating the pattern distally can cancel the descending AP.

**But the distal stimulus fires in both directions too.** Its proximal-going AP collides with the
descending one (both disappear), while its distal-going AP travels on to the muscle. Collision alone
therefore swaps one AP to the muscle for another. The study checks this (§5) and tests the fix the
tutor mentioned, a "reversed anode–cathode configuration" on the distal site:

**Cathodal block.** Around a cathode the activating function has a positive peak (depolarisation,
where the AP starts) flanked by two negative lobes (hyperpolarisation, "virtual anodes"). At a few
times threshold the AP escapes through the lobes in both directions; at a much higher current (the
course slides: above ~8 × threshold) the lobes are hyperpolarised enough to stop it in both
directions. Then the cathode's own AP is trapped, and an AP arriving from the proximal site while the
pulse is on is stopped by the proximal-side lobe. This is the version of the supervisor's scheme that
can cancel the descending AP without sending a new one to the muscle. Its costs: a large current
(safety, and off-target excitation of other fascicles at the distal site) and a timing window (the
pulse must be on when the AP arrives).

**Anodal block.** An anode makes the outside of the membrane positive, which hyperpolarises the
membrane under it (activating function < 0). A strong enough hyperpolarisation keeps the Na channels
from reaching threshold as the AP arrives, so the AP stops there. The anode must be on while the AP
passes, so it needs long pulses (hundreds of µs to ms) and currents several times the excitation
threshold. Next to the hyperpolarised region the membrane is depolarised ("virtual cathodes"), which
can excite the fibre at high currents, and switching the anode off can fire an AP ("anode break").

Other block methods (not in this round): kilohertz high-frequency block (KHFAC), which needs the
frequency-dependent tissue model (`design.md` §5.5).

**Why the FEM alone is not enough.** The FEM gives V and the activating function: where the membrane
is pushed towards firing or away from it. Whether an AP actually starts, travels, collides or is
blocked depends on the membrane dynamics over time, so a fibre model is needed.

## 3. Pipeline

```
results/main3d/UA01_3d.mat            solved 3D model (main_pipeline_3d.m; the viewer shows it)
        │  scripts/export_block_profiles.m  (MATLAB, ~1 min)
        │     src/block/sampleFibreProfiles.m: V along 453 straight fibres, per contact, per 1 A
        │     distal contacts 43–48 = z-mirror of the solved 1–6
        ▼
results/block/UA01_fibre_profiles.h5  + UA01_block_params.json (copy of config/block_params.m)
        │  scripts/run_block_study.py  (Python + NEURON, a few minutes on 7 processes)
        │     src/block/mrg_fibre.py: MRG fibre, extracellular potential = Σ I_c(t) · V_c(z)
        ▼
results/block/UA01_<pattern>_D<d>_*.csv/json/png
```

- **Superposition.** The FEM problem is linear, so a pattern with any currents, polarities and
  timings is a weighted sum of the 1 A single-contact profiles: V(z, t) = Σ_c I_c(t) · V_c(z). One
  "channel" (a stimulator output) is a weighted set of contacts with its own pulse train.
- **Mirror.** Only station −4 is solved. Because the section is the same along z (D-03) and the two
  halves are mirror images, the profile of contact 43 (station +4) is that of contact 1 flipped in z.
- **Fibres.** Straight lines parallel to z on the viewer grid (every 100 µm), inside the fascicles,
  sampled every 50 µm over |z| ≤ 34.5 mm (the domain ends at ±35 mm; end planes insulated since D-17,
  grounded before). Lines through a contact
  are dropped (a fibre cannot cross metal). `results/block/UA01_fibres.png` shows them.
- **MRG fibre.** The double-cable model of McIntyre, Richardson & Grill (2002): nodes of Ranvier with
  fast and persistent Na, slow K and leak; paranodes (MYSA, FLUT) and internodes (STIN) with a
  myelin sheath and a periaxonal space. Rebuilt in Python from MRGaxon.hoc (ModelDB 3810) so that
  the length (68 mm, 59 nodes for D = 10 µm, since D-17; 60 mm and 53 nodes before) and the diameter are parameters. The FEM potential is
  imposed on every compartment as `e_extracellular`.
- **Counting APs.** At the nodes 5 mm from each end of the fibre: an AP at the proximal end = a
  signal to the brain, at the distal end = a signal to the muscle.

## 4. Setup (once per computer)

1. **NEURON.** On Windows `pip install neuron` does not work (PyPI has Linux/macOS wheels only):
   install it with the Windows installer from https://github.com/neuronsimulator/nrn/releases (it
   includes the compiler) and reopen the terminal. On Linux/macOS: `pip install neuron`.
   Check: `python -c "from neuron import h; print(h.nrnversion())"`.
2. **MRG mechanism.** Download ModelDB 3810 (https://github.com/ModelDBRepository/3810) and copy
   `AXNODE.mod` into `results/mrg/model/` (git-ignored), then compile it there:
   `C:\nrn\bin\nrnivmodl.bat` on Windows (gives `nrnmech.dll`), `nrnivmodl` elsewhere.
3. Python packages: numpy, scipy, h5py, matplotlib.

## 5. How to run

```
matlab -batch "setup_paths; export_block_profiles"
python scripts/run_block_study.py --pattern AS2          # patterns are defined in config/block_params.m
python scripts/run_block_study.py --pattern AS2 --diam 16
```

All the numbers (fibre spacing, pulse widths, delays, anodal amplitudes, stimulation level) are in
`config/block_params.m`; re-run the export after changing them, since the Python side reads the JSON
copy.

**Anode–cathode configuration (second round, 2026-11-26): `notebooks/block_anode_cathode.ipynb`.**
Proximal cathode (AS 2 → F7) + distal anode on **the same contacts** (the mirror partners, AS 44:
supervisor 2026-11-26, D-12), with the explanation of each step and of how to read each figure in
the notebook. Free variables: the anodal current and the pulse shape. It maps trapezoidal anodal
pulses (ramps 0–1 ms, against virtual-cathode excitation and anode break) × current on the most
excitable fibre, checks the best settings on all fibres, and runs the textbook bipolar pair
{AS 2 −, AS 44 +}. It uses the same profiles and JSON parameters as the script; the parallel runs
go through `src/block/fibre_jobs.py`. Outputs: `results/block/anode_cathode/`. Run it with *Run
All* (≈ 3–6 min on 7 processes).

**Other distal anode configurations (2026-11-26): `notebooks/block_anode_configurations.ipynb`.**
Our extension beyond the supervisor's rule, after the same-contact anode failed: same proximal
stimulation, and the distal anode chosen from a list (`ANODE_CONFIGS`): the same contact, station +4
(all 6 / 2 cuff pads), the 8 distal cuff pads, the 16 distal shank contacts, all 24 distal contacts (a
"wall"), and the shank contact nearest to the target at each distal station; equal current split.
For each one: a static activating-function screen on all fibres, a ramp × current map on the
example fibre, one raster, the best setting checked on all fibres, and a scoreboard. It needs the
24-contact solve (`p.solve3d.sections = -4:-1`) and a new `export_block_profiles`; with fewer
contacts it skips the configurations it cannot build. Outputs: `results/block/anode_configs/`
(≈ 10–15 min).

**The travelling volley and its anodal block, at fascicle level (third round, 2026-11-27):
`notebooks/block_volley.ipynb`.** The tutor's direction (2026-11-27): the distal site drives the same
contacts as an **anode**; its hyperpolarised region must stop the **global activation** of the target
fascicle, even if some single APs pass. Choices and criterion: D-16. The notebook:
1. computes the threshold of every fibre for the proximal pattern (cached), the recruitment curve of
   every fascicle, and the **maximum selective current** (0.95 × the lowest 10 % threshold of the
   other fascicles); the target is the fascicle the pattern excites most easily (the Gmsh lead
   insertion may have moved the contact to another fascicle);
2. shows the **moving AP** on one fibre (voltage against position and time, snapshots, GIF) and the
   **volley** of each fascicle: the fraction of its fibres carrying an AP at each node and instant
   (`src/block/volley.py`, from the AP times at every node returned by `fibre_jobs.ap_times_job`);
3. sweeps the anodal current × ramp length (trapezoids, plateau 0.1–0.9 ms after the proximal pulse) on
   the target fascicle (proximal + anode, and anode alone to see its own firing), then re-checks the
   best settings on the whole nerve;
4. judges each setting with C1–C3 of D-16 (residual volley at the muscle < 10 % of the target
   fascicle; other fascicles < 10 % at the muscle end and in extra APs to the brain), with the block
   ratio of Tovbis et al. 2026 and the Shannon / CIC check alongside;
5. animates the volley with the block (GIF) and plots the **reach profile** (fraction of the target
   fascicle reached by an AP at each node: it drops past the anode).
It needs the profiles exported from the current (Gmsh) solution first (step B0: `export_block_profiles`;
the notebook refuses older profiles). Outputs: `results/block/volley/` (≈ 20–30 min on 7 processes).
The runs can be replayed in the 3D viewer (*Volley* section, README): `volley.export_viewer_volley` writes the AP
times of each run to `results/viewer/<id>/volley_<name>.js` (a few hundred kB) and the index `volleys.js`.

## 6. Results (first round, 2026-11-25)

Pattern AS 2 alone (intraneural, inside F7, station −4) → distal mirror AS 44 (station +4), 20.6 mm
apart; MRG 10 µm fibres; files `results/block/UA01_AS2_D10_*` (log: `run_AS2_D10.log`).

**Checks.** No AP without stimulus; conduction velocity 51.6 m/s (plausible for a 10 µm MRG fibre);
the distal threshold of each fibre equals the proximal one within 0.98–1.04 (the mirror works; the
small spread is the different node position at the two stations).

**Thresholds (proximal AS 2, 50 µs, 10%-of-fibres criterion).** F7, the fascicle that contains the
contact: 8.6 µA = 0.43 nC. All the other fascicles: 108–225 µA = 5.4–11 nC, the range of the human
TIME thresholds that Zelechowski 2020 reproduced (5.4–6.5 nC). AS 2 is therefore very selective for
F7 (window ×12 before the next fascicle, F8). The F7 value is far below the human data because the
contact sits inside a small fascicle with fibres 100 µm away, with a node aligned on the contact (best
case). Stimulation at I_stim = 1.2 × 8.6 = 10.3 µA excites 5 fibres, all in F7.

**Cathodal block.** On the fibre closest to AS 2, 20–200 µA launches an AP in both directions, but
400 µA launches none: the hyperpolarised regions on either side of the cathode block both branches
(course slides 03: block above ~8 × threshold).

**Collision (same pattern repeated distally, delay 0–0.8 ms): it cannot protect the muscle.** For
every delay, every excited fibre sends exactly one AP to the brain and one to the muscle. The raster
(`UA01_AS2_D10_raster.png`, middle panel) shows why: the distal pulse fires both ways; its upward AP
meets the descending one (they annihilate, so the brain receives only the proximal AP: collision
works), while its downward AP travels to the muscle. For delays longer than the travel time
(20.6 mm / 51.6 m/s ≈ 0.40 ms) the descending AP has already passed and the distal site is
refractory, so again one AP reaches the muscle.

**Anodal block with the same single contact reversed (1 ms): it fails, because the anode excites
first.** With the proximal pulse plus a distal anodal pulse of 10–330 µA (×1–×32), one AP still
reaches the muscle. Run without the proximal pulse, the anode alone already fires every excited fibre
from 21 µA (×2) up: the **virtual cathodes** next to the anode depolarise the fibre at the pulse onset
(right-hand raster: bright regions just above and below the distal site, black = hyperpolarised under
the anode). The AP reaching the muscle is the anode's own. At 660 µA (×64) the brain receives two
APs. No block was found up to 660 µA. For fibres this close to a 200 × 200 µm contact, anodal
excitation needs much less current than block, the known limitation of a single small anode next to
the fibre.

**Cathodal block (the supervisor's scheme with a strong distal cathode): it works on single fibres,
but only at a current that excites the whole nerve.** Grid: distal AS 44 cathodic at 2–64 × I_stim
(21–660 µA), onset 0–0.6 ms after the proximal pulse, 50/200/500 µs (`UA01_AS2_D10_cathodal.png`).
Of the 168 settings, only two stop any of the 5 descending APs:
- 660 µA, 50 µs, onset 0.5 ms blocks 3 of 5. The pulse lands just as the AP reaches the distal site
  (arrival ≈ 0.40 ms + the pulse); third raster panel: the AP stops at the distal site and nothing
  reaches the muscle;
- 83 µA, 500 µs, onset 0.5 ms blocks 1 of 5.

Re-checked on all 453 fibres, the 660 µA setting excites **every fibre of every other fascicle**
(their thresholds are 108–225 µA): 98% of the fibres send an AP to the muscle and an extra one to the
brain. The reason is a ratio: on this contact the block needs ≈ 60 × the F7 threshold, while the
other fascicles fire at only 12–26 × it. A single contact cannot block its own fascicle without
activating the rest of the nerve. (In the heatmap, "success" also counts the F7 fibres the proximal
pulse does not excite and the distal pulse leaves quiet, so it reads higher than the block of the 5
excited fibres.)

**Run time.** Exact thresholds only for the target fascicle: 32 s instead of ~15 min; the whole
study (≈ 4400 fibre simulations, mostly the cathodal grid) takes ≈ 9 min on 7 processes.

**Safety (rough).** 21 µA × 1 ms = 21 nC on a 0.04 mm² contact = 52 µC/cm², already at the
charge-injection limit of Pt; 660 µA × 1 ms = 1.65 mC/cm², far beyond any electrode material. Even if
a block were found at high current, a single small contact could not deliver it safely.

**What this suggests for the next round** (to discuss with the tutor):
0. ask how the supervisor's device avoids the distal AP reaching the muscle: with the same electrodes
   as cathodes, the model gives either collision (the distal AP still reaches the muscle) or
   cathodal block (only at a current that excites the whole nerve). Possible missing pieces: a
   multi-contact distal pattern that focuses the block, a different pulse shape, or a target
   fascicle where the ratio block/off-target is favourable (to test with AS 1 → F2, AS 4 → F1);
1. a **unidirectional distal site**: the distal stimulus must not send its own AP to the muscle, so
   the distal cathode needs an anode on its distal side (e.g. cathode on station +3, anode on +4;
   needs station −3 solved and mirrored). The same virtual-cathode problem must be checked;
2. **anodes spread over several contacts or on the cuff pads** (larger distance from the fibres,
   lower charge density) and **ramped pulses** (a slow rise lets the virtual cathodes accommodate; a
   slow fall avoids anode break, "quasi-trapezoidal" pulses);
3. the fibre-diameter sweep (block and excitation scale differently with D);
4. kilohertz block (later, `design.md` §5.5).

**Gmsh model (D-14, 2026-11-26).** With the new default mesher every station has its own random lead
insertion, so the 48 contacts are all solved and no mirror is used (the export mirrors only unsolved
contacts); fibre lines are placed inside the irregular outlines, and lines that cross a displaced
fascicle edge or the PDMS lead are excluded. The first-round results above are from the PDE Toolbox
model (straight shanks, circles) and must be re-run on the new solution.

**Third round, first run (2026-11-27, Gmsh model, grounded end planes): the anode result is an artefact.**
Files: `results/block/volley_groundEnds/` (kept for the comparison; D-16, D-17).
- *Proximal stimulus.* In the Gmsh model AS 2 lies in the epineurium (lead insertion), 0.14 mm from F7,
  0.16 mm from F3, 0.39 mm from F8: the selectivity window for F7 falls from ×12 (PDE Toolbox model, contact
  inside F7) to **×1.4**. F7 10 % threshold 35.6 µA (1.8 nC); maximum selective current 46.9 µA (2.3 nC)
  recruits 59 % of F7 (13 of 22 fibres) and 7 % of F3. Conduction velocity 51.6 m/s; the volley reaches
  the distal site 0.50 ms after the pulse. The moving AP and the volley figures/GIFs behave as expected.
- *Distal partner.* The "same electrode" AS 44 is not in the same place: each station has its own random
  lead insertion, so AS 44 sits at (0.32, −0.34) mm vs AS 2 at (0.29, +0.24) mm; F7 lies between them
  (0.14 mm from both).
- *Anode sweep (36 settings, 10 µA–2 mA, ramps 0/0.3/1 ms).* No setting met C1; **block ratio 0 at every
  setting**. Up to 16 µA the anode changes nothing; from 26 µA the anode alone fires **all** F7 fibres,
  identically for every ramp. The membrane voltage shows why: the AP starts at the **distal fibre end**
  (z ≈ +30 mm, 5 mm from the grounded end plane), ~0.7 ms after the anode onset, in the nearest and the
  farthest fibre alike, while the anode's own hyperpolarised zone is clearly visible at +10 mm. The anode
  current flows along the nerve to the grounded plane: every fibre sees the same end depolarisation
  (+39.5 V/m² per µA at the end node). Real fibres do not end there, so this firing is an artefact; it hid
  any block above 26 µA, and it also explains the ramp-independent anodal firing of the second round.
- *Fix (D-17):* insulated end planes, fibres to ±34 mm, cache fingerprint with the boundary conditions;
  re-solve, re-export and re-run. The notebook now flags any AP that starts at a fibre end.

**Third round, second run (2026-11-27, insulated end planes, D-17): no fibre-end firing; the anode now fires
the fibres at the end of the insulating sleeve.** Files: `results/block/volley/AS2_monopolar_AS44/`.
- *Near field unchanged*: proximal thresholds, recruitment, conduction velocity and arrival time identical to
  the grounded run (to 3 digits): the end condition does not affect the stimulation side (D-04 comparison).
- *Anode sweep*: still no setting meets C1 and the block ratio is 0 everywhere. Up to 26 µA the anode does
  nothing; from 42 µA (below the proximal current, 46.9 µA: the 0.8 ms anodal pulse needs much less
  current than the 50 µs proximal one) the anode alone fires all of F7, for every ramp. No AP starts at a
  fibre end any more (`frac_fired_at_fibre_end` = 0).
- *Where the anode fires*: the first AP starts at **z = 25.3–26.5 mm**, the end of the 50 mm PDMS sleeve
  (±25 mm), not at the virtual cathodes beside the anode (8.1 and 12.7 mm). Threshold 42 µA: start at
  0.84 ms (late in the plateau); 69 µA: 0.49 ms; 111 µA: 0.34 ms (strength–duration); with ramps it starts
  during the ramp. Inside the insulating sleeve the anode current flows along the nerve and leaves into
  the muscle at the sleeve end: the potential bends there and the activating function has a broad (~5 mm),
  unopposed depolarising peak at 25.4 mm (median 34 V/m² per µA on F7, 40 on F1, 56 on F4: every fascicle).
  Weaker than the virtual cathodes beside the anode (91–302) but not flanked by a hyperpolarised zone, so a
  long pulse fires it first. This is the cuff-edge effect that Tovbis et al. 2026 control with the distal
  electrode edge length (≤ 1 mm in their best design; here the anode is 14.7 mm from the sleeve end).
- *Real or not*: it depends on the device: a tight insulating sleeve 50 mm long is a group assumption (D-05,
  wall and fit unknown). To ask the tutor. A balanced anode (current returned by contacts inside the sleeve
  instead of the far ground) would send almost no net current to the sleeve end.
- *Is a block hidden under the sleeve-edge firing? No.* Proximal APs counted 3.5 mm beyond the anode
  (z = 13.8 mm) at their expected arrival (±0.1 ms): 11/13 at 10–26 µA (the reference: the window misses two
  late fibres; stage 1 shows all 13 reaching the muscle), 9/13 at 42 µA (the sleeve end fires only at
  0.84 ms, after the volley has passed: a clean test), 12–13/13 at 69 µA rectangular (the anode shifts the
  arrival times rather than stopping the APs). The 0/13 cases (111 µA with 0.3 ms ramps, 69–111 µA with 1 ms
  ramps) are collisions with the anode's own AP, fired before the volley arrives, not blocks. **Conclusion:
  below its own excitation threshold the single anode on AS 44 blocks at most 2 of the 13 fired fibres
  (≤ 15 %); C1 needs ≥ 11 of 13.** A single small contact cannot hyperpolarise F7 enough before it (or the
  sleeve end) fires the fibres.
- *Why the timing is not the problem, and what is.* The anode plateau (0.1–0.9 ms after the proximal onset)
  contains the arrival of the volley at the anode (0.49–0.51 ms): the hyperpolarisation is fully on when the
  APs pass. The limit is spatial: a 0.2 mm contact hyperpolarises a region shorter than one internode
  (1.15 mm for 10 µm fibres). On the median F7 fibre the activating function per µA at the nodes is
  +118 (8.1 mm), −90 (9.3), **−655 (10.4, under the anode)**, +7 (11.6), +112 (12.7): only one node is
  strongly hyperpolarised (15 of 22 fibres; two nodes in the other 7). A myelinated fibre conducts with a
  safety margin (the current from an active node is several times what the next node needs), so an AP can
  jump over a single node that cannot fire; it is only delayed (seen at 69 µA, diagnostic 3). To silence
  that one node the current must rise until the virtual cathodes (ratio ~1:4 to the hyperpolarisation) or
  the sleeve end fire the fibres first. Anodes on consecutive stations along F7 (AS 38 at +3, AS 44 at +4)
  widen the hyperpolarised zone to 3 nodes (−196, −72, −245 per µA of total current) but weaken it per µA:
  a trade-off only the fibre model can settle.
- *Next (D-18): balanced anode.* The anodal current is returned by the other 23 contacts of the distal site
  instead of the far ground: the static check predicts the sleeve-end depolarisation to vanish (107 → 9.7
  V/m² per µA) with the F7 hyperpolarisation kept (−655 → −631); the 2 station pads alone would instead push
  F5 hard (619). Notebook `ANODE_RETURN = "distal_site"`, outputs `results/block/volley/AS2_balanced_*`.

**Third round, balanced anode (2026-11-27, D-18): worse.** Files: `results/block/volley/AS2_balanced_AS44_ret-distal-site/`.
Anode AS 44 with its current returned by the other 23 distal contacts. The sleeve-end depolarisation is gone
(9.7 vs 107 V/m² per µA) and no AP starts at a fibre end, but the anode now fires F7 from **26 µA** (32 % of
the fascicle with a rectangular pulse, 82–100 % with ramps; 42 µA for the monopolar anode), and the block
ratio is still **0** at every current. Cause (static AF on F7): the return contacts closest to F7 act as weak
cathodes next to it, AS 38 (station +3, z = 8.3 mm, 0.37 mm from F7) on the node at 8.1 mm where the anode
already has its proximal virtual cathode (+118 → **+171** per µA), and AS 26/33 (stations +1/+2) on the nodes
at 3.5–7 mm (+37 … +104). The hyperpolarisation under the anode is unchanged (−631 vs −655), so the ratio
hyperpolarisation : strongest depolarisation falls from 5.5 : 1 to 3.7 : 1 and the fibres fire earlier.
**Summary of the third round:** monopolar (grounded ends: artefact), monopolar (insulated ends: sleeve-end
firing from 42 µA, ≤ 2 of 13 fibres blocked before it) and balanced (firing from 26 µA, none blocked): in no
configuration does the single mirror anode stop the F7 volley before it fires the fibres itself. Two
mechanisms set the limit: a 0.2 mm contact hyperpolarises one node only (the AP jumps it), and every way of
raising its current makes some other spot depolarise first (virtual cathodes, sleeve ends, return contacts).

**Fourth round (2026-11-27, D-19): the distal configuration designed by optimisation over all contacts.**
`notebooks/block_anode_design.ipynb` (code `src/block/anode_design.py`): for N = 1, 2, 3, 4, 6 active contacts
and K = 1–3 consecutive block nodes, the configuration that maximises block depth / worst depolarisation
(mixed-integer LP on the unit activating functions, D-19), then a NEURON check of the best design per N
(current sweep with the proximal pulse on F7, whole-nerve check with C0–C3). The best block site is between
the two sites (z ≈ 5–8 mm), not under station +4. A chosen design can be replayed with all the volley
pictures in `block_volley.ipynb` (`DIST_PATTERN = "design:<name>"`). Outputs: `results/block/anode_design/`.

*Result (NEURON, 11 designs, rectangular 0.8 ms pulses, 5 µA–1 mA):* **no design blocked a single F7
fibre.** Each one went from 0 % blocked to "the design alone fires all 22 F7 fibres" within one current
step (×1.5): at 17 µA for the 3-contact single-node design, 25–58 µA for the 2- and 3-node designs.
Static numbers at the last current that fired nothing: hyperpolarisation of the block nodes up to
1.5·10⁴ V/m², about half the activating function the 50 µs proximal pulse needs to fire F7 (≈ 3·10⁴ V/m²).
The worst depolarisation anywhere was only 3–8·10³ V/m², which is enough to fire under a 0.8 ms pulse. So
the block needs more current than the current at which the distal site fires the fibres, which is what
Bhadra & Kilgore 2004 report for hyperpolarising DC block (see the fifth round).

**Fifth round (2026-11-27, D-20): is the fibre model wrong? Textbook check and slow-ramp block.**
`notebooks/block_model_check.ipynb`. Literature: Bhadra & Kilgore 2004 (*IEEE TNSRE* 12:313,
doi:10.1109/TNSRE.2004.834205), DC block in NEURON (squid, frog and mammalian models) and in the frog
nerve. Both polarities block, depolarisation block needs less current than hyperpolarisation block, the
block can sit under the electrode or at the virtual electrodes, and a fast switch-on fires the fibre.
Every pulse used here so far rose in ≤ 1 ms. The notebook does three things:
- **A:** reproduces their test on our MRG fibre (point electrode, homogeneous medium, abrupt vs 5 ms ramp,
  anode and cathode). If our fibre blocks there, the model and the code are fine.
- **B:** finds where and when the round-4 designs fired F7: at switch-on, during the plateau (virtual
  cathode) or after switch-off (anode break).
- **C:** tries a slow-ramp anodal block in the nerve (ramps of 2 and 8 ms, plateau while the volley passes),
  with C0–C3 and the charge check. The cost is a long phase, and the charge per phase rises with it.

Outputs: `results/block/model_check/`.

*Results (2026-11-27):*
- **A (point electrode, 10 µm MRG fibre):** our fibre blocks with both polarities, so the model and the
  code can block. With an abrupt switch-on the anode blocks the launched AP from 2.3× the cathodic 50 µs
  threshold, but fires the fibre itself from 1.2–1.7×. A 5 ms ramp does not prevent the firing: it still
  happens during the ramp or the plateau. The ramped cathode blocks from 4.6× and fires from 0.6×. No setting
  gives a clean block. This matches the literature on MRG: Miles et al. 2007 found ramps do not remove the
  onset response. Peña et al. 2025 reproduced the onset-free ramped block only after adding slow sodium
  inactivation, which MRG lacks (`docs/literature/block_methods_review.md`).
- **B (where the round-4 designs fired F7):** never near the block nodes. Every design fired all of F7 at
  the **sleeve end** (z = 25.3–27.6 mm) or at the last fibre node (33.4 mm), either during the plateau or at
  switch-off. Every design had a net anodic current (+0.7 to +2.8 × the largest weight); that current leaves
  the insulating sleeve through its ends. Static check: the sleeve-end depolarisation on F7 is
  ≈ 31 V/m² per µA of net current, spread over 4 nodes. It fired at 1.2–1.4·10³ V/m² in every design, a
  consistent threshold about 20× below the focal 50 µs threshold. The optimiser rated it as harmless (its
  peak is 3–8× lower than the centre peak): **the activating function misses broad depolarisations under
  long pulses.**
- **C (slow-ramp anode, 2/8 ms, single contacts and two designs):** all fire all of F7 from 33–56 µA, with
  no block before. The charge per phase is over the Shannon limit from 55 µA (2 ms ramp).
- **Balanced designs (static check after the run):** with zero net current the sleeve-end bump drops to
  0–18 % of the centre peak. The returns go to proximal contacts (AS 3, 8, 13, 20), and the static ratio
  falls to 0.4–1.3 for K = 2, against 1.4–2.8 unbalanced.

## 7. Limitations and next steps

Final state (review of 2026-11-28). Conclusion of the five rounds of §6: **with this electrode and the
MRG fibre model, no distal configuration tried blocks the volley of the target fascicle before the
distal site fires the fibres itself.** That covers the mirror anode (monopolar and balanced),
optimised multi-contact designs and slow ramps. The fibre model can block (textbook check, round 5 A).
What prevents it is that a block needs more current than the current at which some other point
depolarises (virtual cathodes, the sleeve ends under a net anodic current, return contacts).

Limitations of the study:
- One fibre diameter per run (10 µm); fibres straight and parallel; the same fascicles at both sites (no
  fascicle exchange along the arm, which a real nerve has).
- All 48 contacts are solved, so the first-round mirror trick (section −k → +k) is no longer used. The
  distal profiles come from the solved distal contacts, and the bent leads of the Gmsh model make the
  two sites different anyway.
- Multi-contact patterns use superposition with insulating idle contacts (D-07, approximate).
- Node alignment: a node sits exactly at the proximal station (best case); at the distal sites the
  node position is whatever follows from the internode length.
- Quasi-static tissue, no electrode–tissue interface; the MRG membrane lacks slow sodium inactivation,
  which Peña et al. 2025 needed to reproduce onset-free ramped block (round 5).
- The activating function misses broad depolarisations under long pulses (round 5 B), so the static
  screens and the optimiser (D-19) must always be confirmed in NEURON.
- Safety: the notebooks check the charge per phase and the charge density against the Shannon
  criterion (k = 1.75, `volley.shannon_ok`). The CIC of the contact material (slides 08) is not checked.
  Shannon was fitted on macro-electrodes and is only indicative for 200 × 200 µm contacts.

Options left open (for the tutor / a continuation): balanced designs that also account for the broad
sleeve-end depolarisation; an MRG variant with slow Na inactivation; kilohertz-frequency (KHFAC) block;
fibre-diameter sweep (conduction velocity ∝ D, so the timing window changes). Literature survey of the
block methods: `docs/literature/block_methods_review.md`.

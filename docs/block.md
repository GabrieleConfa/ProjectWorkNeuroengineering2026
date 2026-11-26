# Blocking the signal towards the amputated limb — conduction-block study

Status: first round (2026-11-25), on the solved 3D model UA01 as shown in the viewer. Choices and
their reasons: `decisions.md` D-12. Roadmap item: `design.md` §5.3.

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
  sampled every 50 µm over |z| ≤ 33 mm (the domain ends at ±35 mm, grounded). Lines through a contact
  are dropped (a fibre cannot cross metal). `results/block/UA01_fibres.png` shows them.
- **MRG fibre.** The double-cable model of McIntyre, Richardson & Grill (2002): nodes of Ranvier with
  fast and persistent Na, slow K and leak; paranodes (MYSA, FLUT) and internodes (STIN) with a
  myelin sheath and a periaxonal space. Rebuilt in Python from MRGaxon.hoc (ModelDB 3810) so that
  the length (60 mm, 53 nodes for D = 10 µm) and the diameter are parameters. The FEM potential is
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

## 7. Limitations and next steps

- One fibre diameter per run; fibres straight and parallel; the same fascicles at both sites (no
  fascicle exchange along the arm, which a real nerve has).
- Mirror valid only while the section is the same along z; with Marco's new fascicle shapes or a bent
  shank this must be re-checked (or the distal contacts solved).
- Multi-contact patterns use superposition with insulating idle contacts (D-07, approximate).
- Node alignment: a node sits exactly at the proximal station (best case); at the distal station the
  node position is different, so the "same pattern" does not excite exactly the same fibres.
- Quasi-static tissue, no electrode–tissue interface.
- **Safety not checked yet**: anodal block needs large currents on 200 × 200 µm contacts; the charge
  per phase and the charge density must be compared with the Shannon limit and the CIC of the contact
  material (slides 08) before any of these configurations is considered usable.
- Next: a unidirectional distal site (cathode on one station, anode on a more distal one) needs a
  second solved station (e.g. section −3 → +3 by mirror); multi-contact patterns; fibre-diameter
  sweep; delay window vs fibre diameter (conduction velocity ∝ D).

# Electrical conduction block: literature survey for the volley block (2026-11-27)

Why this note: none of our anodal configurations blocked the F7 volley without first firing the fibres
(`docs/block.md` §6, rounds 3–4). This note collects what the literature says about the ways to block
conduction electrically, and what each means for our intraneural electrode. Sources were found through
PubMed, and only the abstracts were read; none of these papers is in `Papers/` yet. Download a paper and
take numbers from its full text before quoting it in the report.

## 1. Anodal (hyperpolarising) block with quasi-trapezoidal pulses: the tutor's scheme

- **Fang & Mortimer 1991**, *IEEE TBME* 38:168, doi:10.1109/10.76383. Cat medial gastrocnemius nerve,
  **tripolar cuff**, quasi-trapezoidal pulses: square leading edge, **350 µs plateau, exponential trailing
  phase**. The exponential tail **prevents anode-break excitation**. The tripolar cuff **restricts the current
  to the inside of the cuff and so removes the stimulation by the virtual cathode**. Larger fibres are
  blocked at lower current than smaller ones (the blocking threshold correlates with fibre diameter).
- **van den Honert & Mortimer 1979** and **Ungar, Mortimer & Sweeney 1986**: see `tovbis_2026.md`
  (unidirectional APs; plateau ≈ 350 µs, fall ≈ 600 µs).
- **Grünewald et al. 1998**, *World J Urol* 16:329, doi:10.1007/s003450050076. Tripolar cuff and
  quasi-trapezoidal pulses on the sacral roots of dogs: they block the large fibres (sphincter) and spare the
  small ones (bladder). This is "anodal arrest" in a clinical-type application.
- **Bugbee et al. 2001**, *Med Eng Phys* 23:29, doi:10.1016/s1350-4533(01)00013-3. An implantable
  stimulator built for quasi-trapezoidal anodal blocking through nerve cuffs (20 µA–5 mA, up to 1 ms).
- **Vučković, Struijk & Rijkhoff 2005**, *Med Biol Eng Comput* 43:365, doi:10.1007/BF02345814. Simulation
  with a tripolar cuff and long pulses (> 400 µs). With undulating fibres and variable internode lengths,
  the blocking threshold rises by **up to 4×**; it varies by ±60 % (the excitation threshold by ±40 %) and
  differs between the proximal and distal directions. A straight-fibre model (ours) is optimistic.
- **Tovbis et al. 2026**: see `tovbis_2026.md`. FEM + MRG, asymmetric tripolar cuffs, up to 91 % block, at
  currents above the Shannon limit.
- **Ueno et al. 2011**, *IEEE EMBC* 912, doi:10.1109/IEMBS.2011.6090204. Myelinated axon model,
  staircase-shaped stimuli: four response types (subthreshold, cathodic excitation, anodal block, anodal
  break excitation), mapped against the amplitude and the electrode–axon distance.

**For us.** Every anodal-block success in this list uses a **cuff around the nerve**: the insulation
removes the virtual cathode, and an exponential tail removes anode break. Our anode is a 0.2 mm contact
inside the nerve. Its virtual cathodes sit in the same tissue as the fibres, and the 50 mm sleeve has open
ends that act as a virtual cathode (round 3: firing at z ≈ 25 mm). This is the mechanism by which our runs
fail.

## 2. Direct-current (DC) block, ramped DC, charge-balanced DC

- **Bhadra & Kilgore 2004**, *IEEE TNSRE* 12:313, doi:10.1109/TNSRE.2004.834205. NEURON (squid, frog and
  mammalian models) and frog in vivo: **depolarisation block needs less current than hyperpolarisation
  block**; the block can occur under the electrode or at the virtual electrodes; the authors conclude that
  electrical block is generally a depolarisation block, through sodium-channel inactivation.
- **Ackermann et al. 2011**, *Med Biol Eng Comput* 49:241, doi:10.1007/s11517-010-0679-x. Rat sciatic nerve:
  cathodic DC blocked the whole nerve at a mean of 1.73 mA. **Ramped DC blocked without activating the
  nerve**, but down-ramps were more reliable than up-ramps, and the activity depended non-monotonically on
  the up-ramp time. DC also caused a **prolonged conduction failure, likely nerve damage**.
- **Vrabec et al. 2016**, *Med Biol Eng Comput* 54:191, doi:10.1007/s11517-015-1385-5. High-capacitance
  electrodes (platinum black) and a charge-balanced waveform: repeated 13 s DC block of the rat sciatic
  nerve, threshold −0.15 to −3.0 mA; 2 of 6 electrodes still degraded the nerve.
- **Vrabec et al. 2017**, *IEEE TNSRE* 25:517, doi:10.1109/TNSRE.2016.2589541. Charge-balanced DC
  ("CBDC") cycled over 4 contacts ("carousel") gives a continuous block (> 95 %), 22 min without
  degradation.
- **Bhadra et al. 2018** (review), *Bioelectron Med* 1:39, doi:10.2217/bem-2017-0004.

**For us.** Part C of `block_model_check.ipynb` (slowly ramped anode) is supported by Ackermann 2011, which
blocked with ramped DC without activating the nerve. Two caveats come with it. First, safety: DC needs
special high-capacitance electrodes and still risks damage, which our 0.2 mm Pt contacts are not. Second,
the down-ramp matters as much as the up-ramp.

## 3. Kilohertz-frequency AC (KHFAC) block: the established reversible block

- **Bhadra et al. 2007**, *J Comput Neurosci* 22:313, doi:10.1007/s10827-006-0015-5. **MRG model** (axons
  7.3–16 µm), sine waves 1–40 kHz: block between 3 and 40 kHz. The threshold rises with frequency and falls
  with diameter. There is an onset firing at switch-on, and the mechanism is a steady depolarisation of
  several nodes.
- **Kilgore & Bhadra 2014** (review), *Neuromodulation* 17:242, doi:10.1111/ner.12100.
- **Intrafascicular KHFAC: the closest to our device.**
  - **Ackermann et al. 2010**, *Muscle Nerve* 41:117, doi:10.1002/mus.21496: KHFAC block through an
    intrafascicular electrode is feasible (rat).
  - **Dowden, Wark & Normann 2010**, *Muscle Nerve* 42:339, doi:10.1002/mus.21678: Utah slanted arrays in the
    cat sciatic nerve; single electrodes selectively abolished the activity of targeted muscles (conduction
    block at 16 kHz in one animal).
  - **Duncan et al. 2019**, *Neuromodulation* 22:597, doi:10.1111/ner.12802: Utah slanted array in the
    **human median nerve** (one patient with a brachial plexus injury). KHFAC reduced pain sensitivity in one
    finger area while touch and movement were largely preserved, i.e. a **selective intraneural block in a
    human median nerve**.
  - **Regnacq et al. 2025**, *J Neural Eng* 22, doi:10.1088/1741-2552/adc62a: longitudinal intrafascicular
    electrodes (LIFEs), rat: **fascicle-selective, progressive** kHz block, tuned by amplitude and frequency.
    The onset spillover was small and decreased at higher frequency.
- **The onset response and how to remove it.**
  - Miles et al. 2007, *J Neural Eng* 4:390, doi:10.1088/1741-2560/4/4/005: in the MRG model and in rats, a
    linear amplitude ramp from 0 does **not** remove the onset.
  - Vrabec et al. 2019, *JNER* 16:80, doi:10.1186/s12984-019-0554-4: a ramp starting from ≥ 50–75 % of the
    block amplitude avoids a second onset.
  - **Peña, Pelot & Grill 2025**, *J Neural Eng* 22, doi:10.1088/1741-2552/add20e: slow ramps at ≥ 20.8 kHz
    block **without onset** in rats. In axon models the ramp effect appeared only after **adding a slow
    sodium-inactivation gate**, which MRG does not have.
  - Yi & Grill 2020, *PLoS Comput Biol* 16:e1007766, doi:10.1371/journal.pcbi.1007766: a 10 kHz waveform
    optimised to drive the sodium channels into closed-state inactivation blocks without onset (in models).
  - DC + KHFAC: Franke et al. 2014, *J Neural Eng* 11:056012, doi:10.1088/1741-2560/11/5/056012, and Vrabec
    et al. 2013, doi:10.1109/ner.2013.6695927, combine the two with no onset.

**For us.** KHFAC is the method with intraneural, fascicle-selective and human-median-nerve evidence. Its
cost for our model: at ~10–20 kHz the quasi-static purely resistive tissue model is questionable (Horn
2023: the perineurium corner frequency is ≈ 2.6 kHz). The FEM would need a frequency-dependent
(capacitive) perineurium, or at least a check of how much the potential changes. MRG is accepted for KHFAC
(Bhadra 2007).

## 4. What this means for the model ("is NEURON wrong?")

- MRG is the standard model in this literature (Bhadra 2007, Miles 2007, Tovbis 2026). Nothing here
  suggests that it cannot block.
- One known limitation matters for slow ramps: **MRG has no slow sodium inactivation**. Peña et al. 2025
  needed one to reproduce the ramp effects. If our slow-ramp runs still fire, that may be this missing
  mechanism rather than the physics.
- Straight fibres with aligned nodes are the best case. Real fibres undulate, and the blocking threshold
  can be up to 4× higher (Vučković 2005).

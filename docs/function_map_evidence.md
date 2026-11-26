# Which fascicle carries what at mid-arm: evidence (searches of 25 Nov 2026)

Short answer: **there is no published quantitative map of the mid-arm median nerve that tells which fascicle goes to which
muscle.** What exists is listed below with its limits. The model therefore stores a *composition* per fascicle and treats
the assignment of virtual patients as a hypothesis.

| Finding | Source | Reliability for us |
|---|---|---|
| Branches leave in the forearm (exit point, % of forearm length, 0 = wrist, 100 = elbow): pronator teres 98.4 +/- 2.5 (upper) and 88.2 +/- 3.9 (lower), interosseous 81.0 +/- 3.0, FDS 45.8 +/- 4.3, palmar cutaneous (only sensory branch) 16.5 +/- 6.6. Palmaris longus and epitrochlear branching probability reaches 50 % at about 110 and 105 %. About 9 branches per nerve, 71 % muscular. The recurrent (thenar) branch leaves outside the forearm. | Delgado-Martinez 2016 (read in full text) | Good, but n = 3 to 5 nerves. At mid-arm (about 157 %, estimate) all these fibres are still in the trunk. |
| About one eighth of the axons are motor (ChAT): 2836 +/- 414 of 20428 (level 120 %), 2333 +/- 179 of 20028 (level 140 %); 14 and 12 % is our arithmetic. Motor axons are in patches, whose position does not depend on fascicle size or position and cannot be found in advance. | same | Good, n = 5 nerves, two upper-arm levels. |
| Upper arm: the destination of a fascicle cannot be predicted a priori, from size or position. Fascicles reach their organisation in a short segment distal to the epicondyles. | same | Good. |
| Anterior interosseous fascicles form a somatotopic group in the upper arm, at the dorsal and lateral (radial) side of the trunk, about 14.6 +/- 5.4 cm above the humeroradial joint. | MR neurography, Neurology 2014 (PMC3963415, abstract and figure legend) | One group only; patients, not healthy nerves. |
| Position of the anterior interosseous group varies between dorsomedial, dorsolateral and purely dorsal in a dissection of 30 upper limbs. | search snippet, not opened | Shows the variability of even this group. |
| The thenar motor branch varies a lot (more than half of hands deviate from the standard). | search snippet, not opened | Variable. |
| Classical internal topography maps (Sunderland 1945, Jabaley 1980, Williams and Jabaley 1986). | not accessible online (paywalled) | Schematic drawings, no coordinates. To be read from the originals. |
| Hand sensory fibres travel in the same trunk (small fascicles, presumed). The article also suggests that intrinsic hand fibres may sit in larger central fascicles: the two statements are not reconciled. | Delgado-Martinez 2016 | Weak. |

**Not found:** positions of the pronator teres, FCR, FDS, thenar and sensory groups at mid-arm; a histogram of fibre diameters
for motor and sensory fibres of the median nerve at mid-arm; sizes of every fascicle of one nerve together with its function.

**Consequence for the model:** the function of a fascicle is stored as a composition (fractions of fibres per group), set by
a data file for a measured patient and drawn from a stated hypothesis for a virtual patient (`assignFunctions`: random,
forearm rule, MR neurography). The first useful measurement is a real section with orientation and, ideally, a map of the
fascicle groups obtained by stimulation or by tracing.

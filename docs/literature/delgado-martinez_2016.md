# Delgado-Martínez et al. 2016 — Fascicular topography of the human median nerve

**Citation.** Delgado-Martínez I, Badia J, Pascual-Font A, Rodríguez-Baeza A, Navarro X.
*Fascicular Topography of the Human Median Nerve for Neuroprosthetic Surgery.* Front. Neurosci.
10:286, 2016. doi:10.3389/fnins.2016.00286. File: `Papers/fnins-10-00286.pdf`.

**Role in the project.** Main source for **phase 1 (geometry)** and for the **virtual-patient
variability** of phase 6: it gives the number, size and position of the fascicles of the human
median nerve at different arm levels, and where an intraneural electrode should be implanted.

## What they did

- 8 median nerves from 5 human donors. Dissection to map where each muscular branch leaves the
  main trunk; distances normalised to forearm length (**nFD**: 0% = styloid processes at the wrist,
  100% = epicondyles at the elbow; > 100% = upper arm).
- One nerve cut into 11 blocks (I = wrist … XI = upper arm), 950 slices of 25 µm every 250 µm,
  fascicles segmented automatically (Otsu) → 18,122 fascicle cross-sections, 3D reconstruction.
- ChAT immunostaining (marks motor axons) at nFD = 120% and 140% (upper arm).

## Numbers useful for the geometry generator

Morphometry (one nerve, mean ± SEM over blocks; p. 6–7, Fig. 5, Table 1):

| Quantity | Value |
|---|---|
| "Epineurium area" (area of the ellipse enclosing all fascicles ≈ nerve cross-section) | 8.41 ± 1.01 mm² overall; upper arm 10.40 ± 1.75; proximal forearm 7.13 ± 0.59; distal forearm 5.80 ± 0.10; block I (wrist) 15.44 |
| Eccentricity "of the nerve section" | 0.65 ± 0.02 per the text — Fig. 5C plots it as *fascicle* eccentricity: open question (issue D1), see Caveats |
| Endoneurium area (sum of fascicle areas) | 3.12 ± 0.24 mm² overall; upper arm 4.13 ± 0.23; forearm 2.46–2.73 |
| Endo/epi area ratio | 0.42 ± 0.03 (0.24 – 0.68) |
| Fascicles per section | 19.85 ± 1.23 overall; forearm 21.87 ± 0.58; upper arm (blocks X, XI) 12.81 ± 0.73 and 11.81 ± 0.32 |

Fascicle clusters (k-means on area and distance from the nerve centroid; Table 1):

| Cluster | Area (mm²) | Share | Where |
|---|---|---|---|
| A (small) | ≈ 0.12–0.13 | 94% | everywhere |
| B (medium) | ≈ 0.52–0.55 | ≈ 5% (1.6% in forearm, 18.7% in upper arm) | mostly upper arm |
| C (large) | ≈ 1.5–2.0 | 0.8% | upper arm only |

Distance classes from the nerve axis: 0 ≈ 0.68 mm, 1 ≈ 1.05–1.46 mm, 2 ≈ 2.1–2.35 mm.
Fascicles are spread homogeneously over the section regardless of size (Fig. 8).

Axons (upper arm, p. 7–8): ≈ 20,000 axons per section, of which ≈ 2,300–2,800 ChAT+ (≈ 1/8
motor). Motor axons are **clustered in patches** inside mixed fascicles, at positions not
predictable from fascicle size or position; purely motor or purely sensory fascicles are rare.

## Anatomy / implantation conclusions

- Muscular branches (≈ 6 per nerve) to PT, FCR, PL, FDS, FDP, FPL, PQ (forearm flexors and
  pronators); interosseous nerve leaves at nFD ≈ 81%. Fascicles stay compact near the axis until
  their exit point; forearm-muscle fascicles lie dorsally, around the hand fascicles.
- Three regions: **proximal** (upper arm; few large fascicles, motor axons mixed and
  unpredictable), **early distal** (antecubital fossa, ≈ 50 mm below the epicondyles; topography
  becomes defined; all motor fibres for the hand still present), **late distal** (forearm to
  carpal tunnel; many small fascicles, mostly hand/sensory).
- **Intrafascicular electrodes (TIME, LIFE)**: best in the upper arm and the proximal third of
  the forearm (the 95% CI of the interosseous branching point, nFD ≈ 77–89%). Cuff/FINE: distal
  forearm.

## How we use it

- Choose the **implantation level** first (upper arm vs proximal forearm): it fixes fascicle
  count (≈ 12 vs ≈ 22) and size mix (B/C clusters vs almost all A). Log the choice in
  `docs/decisions.md`.
- Generate virtual patients by sampling nerve area, eccentricity, fascicle number and areas
  from the statistics above (seeded RNG), placing non-overlapping fascicles inside the ellipse.

## Newer literature (checked 2026-10-21)

No later study replaces this one (3D fascicle tracking + motor axons along the whole arm). Newer
data, summarised in `Papers/main.tex` (Sec. "What has been published since 2016"):
`brill-tyler_2017.md` (4 nerves, outer area, shape ratio 0.51 ± 0.13, counts with SD),
`abdalbary_2021.md` (perineurium thickness, doubtful), `romeni_2023.md` (reuses these sections
for TIME optimisation; grounded 10 mm saline cylinder). Also cited there: Tovbis 2020 (3D
reconstruction of 5 forearm median nerves), Cartwright 2017 (in-vivo 70 MHz ultrasound: 22.68
fascicles at the wrist, 40 wrists), Yoon 2018 (7T MRI resolves fascicles), Upadhye 2022 (vagus
fascicles split/merge every ≈ 560 µm).

## Caveats

- Morphometry comes from **one** nerve; inter-subject variability of the cross-section is not
  quantified (only of branching points, n = 5). Any spread we use for virtual patients is partly
  an assumption — say so in the discussion.
- Values are mean ± **SEM**, not SD: the population spread is larger than the ± suggests.
- Paraffin-embedded fixed tissue shrinks; in-vivo dimensions are likely larger.
- **Nerve eccentricity is not really reported** (seen 2026-10-19): the only eccentricity plot
  (Fig. 5C) is labelled "Fascicle Eccentricity" and its 11 bars average ≈ 0.65, i.e. the number
  the text attributes to the nerve section; the text also says the epineurium eccentricity
  "remained unchanged (data not shown)". Open question for the tutor (issue D1 in
  `Papers/main.tex`, Appendix A); the schematic outlines in Fig. 4 look like ellipses with axis
  ratio ≈ 1.3–2.
- The same paragraph gives the endo/epi ratio block by block (0.24–0.68) and then says it
  "remained unchanged" (issue D2).
- The 95% CI of the interosseous branching point (76.67–89.32% nFD) is not centred on the mean
  81.0 ± 3.0 (SEM, n = 5 → t-interval 72.7–89.3); the two PT intervals are consistent (issue D3).
- The paper says the interosseous nerve also supplies FDS (and FCR in 3/5); textbooks give only
  FPL, lateral FDP, PQ (issue X4).
- Fig. 8C shows distance class 2 spanning ≈ 1.9–2.8 mm, and Fig. 8A distances up to ≈ 4 mm.
- [CHECK] Distance class 2 (≈ 2.2–2.35 mm from the centroid) is larger than the semi-major axis of
  an 8.4 mm², e = 0.65 ellipse (≈ 1.9 mm); probably driven by the larger proximal blocks. Check
  against Fig. 8 before using the distances directly.

# Brill & Tyler 2017 — Quantification of human upper extremity nerves and fascicular anatomy

**Citation.** Brill NA, Tyler DJ. *Quantification of human upper extremity nerves and fascicular
anatomy.* Muscle Nerve 56(3):463–471, 2017. doi:10.1002/mus.25534. Open access: PMC5712902
(not yet in `Papers/`).

**Role in the project.** Second source for **phase 1 (geometry)** and the **virtual patients**
(phase 6): unlike Delgado-Martínez 2016 (one nerve), it gives a between-subject SD and the
*outer* nerve boundary and shape. Summarised in `Papers/main.tex`, Sec. "What has been
published since 2016".

## What they did

- 5 arms, 3 embalmed cadavers (all > 50 y), 4 median nerves (also 4 ulnar, 3 radial); paraffin,
  1-cm blocks; 65 median-nerve cross-sections, 1373 median fascicles.
- Length normalised: 0 axilla, 0.5 elbow (medial epicondyle), 1 wrist (radial styloid);
  9 regions — 2–8 are 1/8 each, 1 and 9 are 1/16. Region 2 = upper arm, 5 = elbow,
  6 = proximal forearm, 7 = mid forearm.
- Nerve area traced on the **outer epineurium**; values are mean ± **SD across nerves**.

## Numbers (median nerve; Tables 2–3, Results)

| Quantity | Value |
|---|---|
| Nerve area (mm²), regions 2 / 5 / 6 / 7 | 15.7 ± 6.78 / 27.59 ± 18.87 / 10.94 ± 2.77 / 8.27 ± 2.18 (n = 3) |
| Minor × major diameter (mm), region 7 | 2.3 ± 0.27 × 4.47 ± 0.71 |
| Minor/major ratio (all sections) | 0.51 ± 0.13 → eccentricity ≈ 0.86 (our computation); flattest at elbow |
| Fascicle count, regions 2 / 4 / 5 / 6 / 7 | 14 ± 2.6 / 18 ± 1.9 / 23.2 ± 3.8 / 19 ± 2.4 / 18.6 ± 3.3; max 42 |
| Fascicle diameter (circle-equivalent) | upper arm (regions 1–5) 0.6 ± 0.3 mm; forearm 0.47 ± 0.18 mm |
| Fascicle circularity | 97% > 0.75 (all nerves) → circles are fine |
| Inner fascicles (not on convex hull) | ≈ 0.58 (N − 6.08) |
| AIN branch point | 3.6 ± 2.1 cm below the elbow (≈ 84% nFD, our computation from Table 1) |

## Caveats

- 3 elderly embalmed donors; no motor/sensory identification; no fascicle-to-target map.
- Shrinkage: they quote 10–43% loss of fascicular area with fixation (Dyck) and a linear
  correction factor 1.25 for paraffin-embedded porcine fascicles (Kundu).
- Text says an exponential fit of mean diameter vs count, the Fig. 5 caption says a power fit
  (minor; not used).

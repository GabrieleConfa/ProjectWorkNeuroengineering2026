# Tovbis, Koh, Jeong, Yoo 2026 — Asymmetric tripolar cuff designs for optimising anodal block

**Citation.** Tovbis D, Koh RGL, Jeong R, Yoo PB. *Asymmetric tripolar nerve cuff electrode designs
for optimizing anodal block of myelinated fibers.* BioMedical Engineering OnLine 25, 2026.
doi:10.1186/s12938-026-01579-6 (PMID 42177531, open access PMC13374333). Not in `Papers/`: read
from the PubMed Central full text on 2026-11-27 (download the PDF into `Papers/` before quoting
page numbers in the report).

**Role in the project.** The closest published work to our block study: **FEM + MRG fibres**
(as we do, D-12) used to find electrode layouts that stop the distal-going APs with an anode, with an
in-vivo check in the rat. It gives us the **block metric** and an idea of the **block level one can
expect** (D-16).

## What they did

- Rat sciatic nerve with a flat-interface nerve electrode (FINE), then a human-size cylindrical cuff;
  MRG fibres 5.7–16 µm (several diameters); proximal cathode, one or two distal anodes (bipolar and
  asymmetric tripolar designs, with a current split between primary and secondary anode).
- For each configuration the current is increased and every fibre is classified: below threshold,
  activated, **blocked** (its distal-going AP stopped at the anode), **reactivated** (the anode's
  virtual cathode fires it again at higher current). Pulse widths 0.5 and 1 ms.
- In vivo (rat): foot-muscle EMG (normalised RMS) as the readout of what passes the block.

## What we take

- **Block metric:** *anodic block* = percentage of the **activated** fibres whose distal-going AP is
  stopped. We report it as `block_ratio` (`src/block/volley.py`, notebook `block_volley.ipynb`).
- **Expected level: partial, never 100 %.** Rat model: maximum block 88 % (computational) and 75 % (in
  vivo); human-size model: up to 91 % (bipolar with a short distal electrode edge, ≤ 1 mm) and 90 %
  (tripolar with a 6:4 current split between the two anodes). The short bipole reached only 18 %: it
  failed to block small and medium fibres (5.7–8.7 µm) far from the anode. This supports the tutor's
  remark (2026-11-27) that some APs pass the block, and our fascicle-level criterion (D-16).
- **Mechanism of failure:** reactivation by the anode's **virtual cathode** at high current; suppressed
  by a secondary anode placed beyond the primary one, or by a short distal electrode edge (the field
  outside the cuff shapes the virtual cathode). Relevant to our single-contact anode, which in the second
  round fired the fibre before blocking it (`docs/block.md` §5).
- **Safety:** they compare the charge and charge density per phase with the **Shannon criterion with
  k = 1.75 for platinum** (citing Shannon 1992) and find that the currents needed for peak block are
  above the limit for the anodes they tested (at 0.5 ms: safe up to 1.6 mA for a 0.1 mm anode,
  3.7 mA for 0.5 mm, 5.0 mA for 1 mm). We use the same k in the rough safety check (`parameters.md`).
- **Node position:** in vivo the nodes of neighbouring fibres are not aligned; they checked the effect
  of the longitudinal node position on the thresholds. Our runs align a node with the proximal contact
  (best case), a limitation to keep.

## Caveats

- Cuff/FINE around the nerve, not intraneural contacts: our anode is a 200 × 200 µm contact inside or
  next to a fascicle, much closer to the fibres; virtual-cathode excitation is stronger there.
- Mostly monofascicular geometry (they discuss it as a limitation); EMG as the in-vivo readout.
- We have not yet read the supplementary material (table T1 on "bidirectional" and "distal" block).

## Related (abstract only, PubMed)

- van den Honert & Mortimer 1979, *Science* 206:1311, doi:10.1126/science.515733 — unidirectional APs
  by anodal block with a quasi-trapezoidal pulse and a tripolar cuff; proposed for collision block of
  motor nerves.
- Ungar, Mortimer & Sweeney 1986, *Ann Biomed Eng* 14:437, doi:10.1007/BF02367364 — quasi-trapezoidal
  pulse minimising the charge: plateau ≈ 350 µs, exponential fall ≈ 600 µs (90–10 %); critical current
  density for arrest 0.47 ± 0.08 mA/mm² in the cuff.
- Sweeney & Mortimer 1986, *IEEE TBME* 33:541, doi:10.1109/TBME.1986.325818 — asymmetric two-electrode
  cuff for unidirectional APs.
- Vučković, Tosato & Struijk 2008, *J Neural Eng* 5:275, doi:10.1088/1741-2560/5/3/002 — pig vagus:
  anodal block suppressed only 0–40 % of the largest-fibre activity, attributed to the large nerve
  (fibres at very different distances from the electrode); charge per phase above the safe limits.
- Ahmed et al. 2020, *Sci Rep* 10:9221, doi:10.1038/s41598-020-66332-y — anodal block makes vagus
  stimulation directional in the rat, judged by functional markers (breathing vs heart rate) and
  compound potentials; both polarities still gave mixed responses (partial block).

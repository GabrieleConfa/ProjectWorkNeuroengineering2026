# Provided electrode layout

`sites.mat` — the 48 active sites (AS) of the new electrode, provided by the tutor (Luca Sbarsi)
together with the description script `src/MatlabCode/load_active_sites.m`. This is **input
data**: do not edit or overwrite it.

## Origin of the file (2026-11-08)

The file arrived inside a macOS zip. Once unzipped, `src/MatlabCode/sites.mat` was a **folder**
containing the real `sites.mat` (491 bytes) and a `__MACOSX/._sites.mat` resource-fork file
(macOS metadata, no data). The real file was copied here unchanged
(MD5 `a0721df9e15b67883225fe0e0607aa50`), the metadata file was discarded, and the empty
folder was removed.

Note: the tutor's `load_active_sites.m` calls `load("active_sites.mat", "sites")`, a file name
that does not exist. The script was left untouched as the original reference; the project code
uses `src/electrode/loadActiveSites.m`, which reads this file, converts it to SI units and
fills in the cuff radius. `tests/testLoadActiveSites.m` checks the layout described below.

## Content

Variable `sites`, 48 rows (one per AS) x 10 columns, **lengths in mm, angles in degrees**:

| Col | Name | Meaning |
|---|---|---|
| 1 | Electrode | AS number, 1...48 |
| 2 | Type | 1 = intraneural (on the shank), 2 = cuff (on the nerve surface) |
| 3 | Section | -4...-1 (z < 0), +1...+4 (z > 0) |
| 4 | z | centre position along the nerve axis |
| 5 | r | distance of the centre from the axis; NaN for cuff sites → use r = R + thickness/2 |
| 6 | theta | angular position of the centre, from x, counterclockwise |
| 7 | phi | direction of the normal to the square face: n = [cos phi, sin phi, 0] |
| 8 | thickness | box size along n |
| 9 | sideZ | face side parallel to the nerve axis |
| 10 | sideXY | face side in the x-y plane, along w = [-sin phi, cos phi, 0] |

Frame: z = nerve axis, origin at the centre of the device; x = r cos(theta), y = r sin(theta).

## What the data say about the device

- 8 sections at z = ±4.3, ±6.3, ±8.3, ±10.3 mm (2 mm apart); section +k is the mirror image of
  section −k with respect to z = 0. The device spans |z| ≤ 10.4 mm (centre ± sideZ/2).
- Each section: 4 intraneural AS on a shank along a diameter, at r = 0.3 and 0.9 mm on both
  sides of the axis (fixed by the device, independent of the nerve radius), plus 2 cuff AS on
  the nerve surface, at 90° from the shank, with the normal pointing to the axis.
- The shank direction rotates by 45° from one section to the next: theta = 0°, 45°, 90°, 135°
  for sections −4, −3, −2, −1 (and +4, +3, +2, +1). So there are only **4 distinct shank
  orientations**, each present twice (at ±z).
- All AS are equal: thickness 0.005 mm, faces 0.2 mm x 0.2 mm.
- Material: gold (tutor's slides, PresentazionePW5_1105.pdf).

## Open points (also in docs/project_brief.md, "Still unknown")

- **[CHECK] Thickness:** the matrix gives 0.005 mm = 5 µm, the tutor's slides give
  "50 µm x 200 µm x 200 µm". One of the two is off by a factor of 10 — ask the tutor.
- **Shank body:** only the AS are given, not the shank substrate that carries the intraneural
  sites (its width, thickness, material). An insulating shank would shield one side of each AS.
- **Which faces are active:** the data give boxes; it is not stated whether current leaves
  through both square faces or only one (cuff sites: presumably only the face on the nerve).

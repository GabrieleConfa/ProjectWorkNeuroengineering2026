# Pelot, Behrend & Grill 2019 — Parameters used in FEM of compound peripheral nerves

**Citation.** Pelot NA, Behrend CE, Grill WM. *On the parameters used in finite element modeling
of compound peripheral nerves.* J. Neural Eng. 16(1):016007, 2019. doi:10.1088/1741-2552/aaeb0c.
Open access: PMC7309635. Not in `Papers/`.

**Role in the project.** Phase 2 (materials) and phase 3 (perineurium representation): explains
where the usual perineurium numbers come from and how temperature/frequency change them; supports
the thin-layer (sheet-resistance) treatment. Summarised in `Papers/main.tex`, Sec. "What has been
published since 2023" (Horn chapter).

## Numbers

Perineurium sheet resistance derived from Weerasuriya 1984 (frog), their Table 2:

| Case | R_s (Ω·m²) | R_s (Ω·cm²) |
|---|---|---|
| DC, room temperature (as published) | 0.0478 | 478 |
| DC, 37 °C (Q10 correction) | 0.0250 | 250 |
| 10 kHz, room temperature | 0.0168 | 168 |
| 10 kHz, 37 °C | 0.0088 | 88 |
| "Common value in literature" (300 Ω·m × 50 µm) | 0.0149 | 149 |

- DC, 37 °C spread over 3% of the frog fascicle diameter → ρ_peri = **1149 Ω·m** (their default).
  Zelechowski 2020 uses 1136 Ω·m; Pelot lists Raspopovic 2017 (Zelechowski's source) among papers
  using the Q10 correction → very probably the same origin (our reading).
- Endoneurium: longitudinal 1.75 Ω·m, transverse 12 Ω·m (Ranck & BeMent 1965) as defaults; bulk
  transverse ≈ 10 × longitudinal by their own 2D estimate.
- Sensitivities: thresholds up to ≈ +130% for ρ_peri ≈ 400 → 2200 Ω·m; +35–250% for
  ρ_endo,long 3.5 → 0.75 Ω·m; +30–180% for larger fascicle radius; nearly insensitive to
  ρ_endo,trans.
- Contact-impedance boundary ≈ meshed perineurium (similar thresholds).

## Caveats

- Cuff electrodes, not intraneural. Derived values rest on frog data and on a Q10 factor.

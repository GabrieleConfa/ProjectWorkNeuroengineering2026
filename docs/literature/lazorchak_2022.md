# Lazorchak et al. 2022 — Capacitive laminae and in-vivo strength–duration thresholds

**Citation.** Lazorchak N, Horn MR, Muzquiz MI, Mintch LM, Yoshida K. *Accurate simulation of cuff
electrode stimulation predicting in-vivo strength-duration thresholds.* Artif. Organs
46(10):2073–2084, 2022. doi:10.1111/aor.14374. Open access: PMC9529855.

**Role in the project.** Evidence for issue Q2 (resistive vs capacitive model). Same group as
Horn 2023.

## What they did

Time-dependent FEM of a bipolar cuff on rat vagus (2 rats in vivo), MRG fibres, pulse widths
30–1000 µs. 8 cases: saline vs in-air outer environment × low/high-frequency resistivity ×
ε_r = 80 vs measured ε_r.

Parameters (their Table 1): epineurium ρ0 = 68, ρ∞ = 1.3 Ω·m, ε_r,measured = 3.7e6; perineurium
ρ0 = 3460, ρ∞ = 273 Ω·m, ε_r,measured = 1440; endoneurium ρ = {12, 12, 1.7} Ω·m; saline 1 S/m.
The measured permittivities come from unpublished ex-vivo canine vagus data.

## Findings

- Best fit r² = 0.96: in-air (as the experiment), low-frequency ρ, measured ε_r.
- Environment and resistivity dominate; permittivity matters mainly at short pulse widths.

## Caveats

n = 2 rats; best-fit model fibre conducts at ≈ 64 m/s vs ≈ 23 m/s for the in-vivo peak;
agreement described as "within a factor of 10". Cuff, not intraneural.

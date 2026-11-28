function safety = checkSafety(thr, sites, params)
%CHECKSAFETY Tissue and electrode safety of the stimulation at the threshold charges.
%
%   NOT IMPLEMENTED (final status, 2026-11-28). Skeleton from the first work plan (docs/work_plan.md),
%   never filled in: the project did this job elsewhere. Calling it stops with an error. Instead:
%   safety is checked by src/safety/fieldExposure.m + exposedVolume.m (field exposure,
%   scripts/run_field_safety.m) and by the Shannon charge check in src/block/volley.py (shannon_ok).
%   Kept only as a record of the planned interface; not used by any pipeline.
%
%   safety = checkSafety(thr, sites, params)
%
%   Owner: WP5. Status: SKELETON.
%
%   Tissue:    Shannon criterion, log10(D) = k - log10(Q), with Q the charge per
%              phase and D = Q / area the charge density (units of the original paper:
%              uC and uC/cm^2 — convert from SI here and only here).
%   Electrode: D below the charge-injection capacity of gold (params.safety.cicGold_Cpm2).
%   With a 0.2 x 0.2 mm face, 6 nC -> 15 uC/cm^2; with both faces active, half of it.
%
%   Inputs
%     thr     table from thresholdPerFascicle (Qth_C)
%     sites   table from loadActiveSites (face sizes)
%     params  struct from default_params() (params.safety, params.electrode.activeFaces)
%
%   Output: table, one row per (site, fascicle)
%     Q_C, chargeDensity_Cpm2, shannonK, tissueOK, electrodeOK

arguments
    thr table
    sites table
    params (1,1) struct
end
error('project:notImplemented', 'checkSafety: not implemented yet (WP5).');
end

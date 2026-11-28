function thr = thresholdPerFascicle(af, fibres, siteId, params)
%THRESHOLDPERFASCICLE Threshold current and charge of each fascicle for one active site.
%
%   NOT IMPLEMENTED (final status, 2026-11-28). Skeleton from the first work plan (docs/work_plan.md),
%   never filled in: the project did this job elsewhere. Calling it stops with an error. Instead:
%   thresholds are computed with the MRG fibre model (src/block/fibre_jobs.py threshold_job) and the
%   10 %-of-fibres threshold per fascicle in notebooks/block_volley.ipynb (section 5.2).
%   Kept only as a record of the planned interface; not used by any pipeline.
%
%   thr = thresholdPerFascicle(af, fibres, siteId, params)
%
%   Owner: WP4. Status: SKELETON.
%
%   By linearity, fibre i fires when I * afMax_i >= afThreshold(D_i), so its
%   threshold current is I_i = afThreshold(D_i) / afMax_i. The fascicle threshold
%   is the current that recruits params.activation.recruitmentFraction of its fibres.
%
%   Inputs
%     af       struct from activatingFunction (for AS siteId)
%     fibres   table from makeFibres
%     siteId   AS number
%     params   struct from default_params() (activation, stim.pulseWidth_s)
%
%   Output: table, one row per fascicle
%     siteId, fascicle
%     Ith_A         threshold current [A]
%     Qth_C         threshold charge per phase = Ith_A * pulse width [C]
%     recruitment   [I_A, fraction] curve (n x 2)
%
%   Validation target: 5.39 +- 0.98 nC, TIME in the human median nerve (Zelechowski 2020).

arguments
    af (1,1) struct
    fibres table
    siteId (1,1) double
    params (1,1) struct
end
error('project:notImplemented', 'thresholdPerFascicle: not implemented yet (WP4).');
end

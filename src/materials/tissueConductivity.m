function c = tissueConductivity(tissueCode, fascicleRadius_m, params)
%TISSUECONDUCTIVITY Conductivity tensor (diagonal) of a tissue, in the createpde "c" format.
%
%   NOT IMPLEMENTED (final status, 2026-11-28). Skeleton from the first work plan (docs/work_plan.md),
%   never filled in: the project did this job elsewhere. Calling it stops with an error. Instead:
%   the conductivities are set in config/pw6_3d_params.m (p.materials) and assigned to the elements in
%   src/main/build3DModel.m (tissueSigma).
%   Kept only as a record of the planned interface; not used by any pipeline.
%
%   c = tissueConductivity(tissueCode, fascicleRadius_m, params)
%
%   Owner: WP5. Status: SKELETON.
%
%   Inputs
%     tissueCode        1 endoneurium, 2 perineurium, 3 epineurium, 4 medium, 5 electrode
%     fascicleRadius_m  endoneurium radius of the fascicle [m] (needed only if the
%                       perineurium is rescaled to keep rho*t constant; NaN otherwise)
%     params            struct from default_params() (params.materials)
%
%   Output
%     c   3x1 [sigmaT; sigmaT; sigmaL] in S/m (x, y across the nerve, z along it),
%         to be passed to specifyCoefficients(model, c = c, Cell = ...)
%
%   Must stop with an error if a required value in params is NaN.

arguments
    tissueCode (1,1) double {mustBeInteger, mustBePositive}
    fascicleRadius_m (1,1) double
    params (1,1) struct
end
error('project:notImplemented', 'tissueConductivity: not implemented yet (WP5).');
end

function sol = solveUnitCurrent(nerveModel, siteId, params)
%SOLVEUNITCURRENT Potential in the nerve when one active site injects 1 A.
%   sol = solveUnitCurrent(nerveModel, siteId, params)
%
%   Owner: WP3. Status: SKELETON.
%
%   Solves div(sigma grad V) = 0 (quasi-static volume conductor) with:
%     - sigma per cell from tissueConductivity (anisotropic: c = [sT; sT; sL]);
%     - current density I/area on the faces of AS siteId (Neumann, "g" term);
%     - outer boundary from params.domain.outerBC (V = 0 or insulated);
%     - all other AS: floating / insulated (to decide and log).
%   The problem is linear: results for any current are I * sol (principle 5).
%
%   Inputs
%     nerveModel  struct from meshNerveModel
%     siteId      AS number, 1...48
%     params      struct from default_params()
%
%   Output: struct (docs/architecture.md §4)
%     siteId, I_A (= 1), result (PDE results object), V_site_V, meshInfo
%
%   Tests to write first: point source in an infinite homogeneous medium
%   (V = I/(4 pi sigma r)), linearity, mirror symmetry between AS of sections +k / -k.

arguments
    nerveModel (1,1) struct
    siteId (1,1) double {mustBeInteger, mustBeInRange(siteId, 1, 48)}
    params (1,1) struct
end
error('project:notImplemented', 'solveUnitCurrent: not implemented yet (WP3).');
end

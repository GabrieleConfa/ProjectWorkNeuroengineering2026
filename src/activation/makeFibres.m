function fibres = makeFibres(patient, params, seed)
%MAKEFIBRES Place myelinated fibres inside each fascicle and compute their nodes of Ranvier.
%   fibres = makeFibres(patient, params, seed)
%
%   Owner: WP4. Status: SKELETON.
%
%   Inputs
%     patient  struct from makePatient (fascicle positions and radii)
%     params   struct from default_params() (params.activation, params.domain.length_m)
%     seed     integer for rng(seed): positions, diameters, node offset along z
%
%   Output: table, one row per fibre
%     fascicle      fascicle index
%     x_m, y_m      position in the cross-section [m] (inside the endoneurium)
%     diameter_m    fibre diameter D [m]
%     internode_m   internodal length = params.activation.internodeFactor * D [m]
%     nodeZ_m       z of its nodes of Ranvier [m] (cell array of row vectors),
%                   with a random offset in [0, internode) so that nodes are not aligned

arguments
    patient (1,1) struct
    params (1,1) struct
    seed (1,1) double
end
error('project:notImplemented', 'makeFibres: not implemented yet (WP4).');
end

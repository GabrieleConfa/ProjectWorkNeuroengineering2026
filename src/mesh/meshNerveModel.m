function nerveModel = meshNerveModel(geom3, patient, sites, params)
%MESHNERVEMODEL Create the PDE model on the 3D geometry and generate the tetrahedral mesh.
%
%   NOT IMPLEMENTED (final status, 2026-11-28). Skeleton from the first work plan (docs/work_plan.md),
%   never filled in: the project did this job elsewhere. Calling it stops with an error. Instead:
%   meshing is done in src/main/build3DModel.m (meshPersonalDevice + mesh_personal.py for Gmsh,
%   buildNerveGeometry3D for the PDE Toolbox mesher).
%   Kept only as a record of the planned interface; not used by any pipeline.
%
%   nerveModel = meshNerveModel(geom3, patient, sites, params)
%
%   Owner: WP2. Status: SKELETON.
%
%   Inputs
%     geom3    struct from addActiveSites
%     patient  struct from makePatient;  sites  table from loadActiveSites
%     params   struct from default_params() (params.mesh)
%
%   Output: struct (docs/architecture.md §4)
%     pde           PDEModel from createpde() with geometry and mesh
%                   (createpde, not femodel: femodel accepts only isotropic sigma)
%     cellTissue, cellFascicle, siteFaces, siteArea_m2   (from geom3)
%     patient, sites                                       (copies of the inputs)
%
%   Mesh sizes per tissue via generateMesh(..., Hcell = {...}), finer at the AS.

arguments
    geom3 (1,1) struct
    patient (1,1) struct
    sites table
    params (1,1) struct
end
error('project:notImplemented', 'meshNerveModel: not implemented yet (WP2).');
end

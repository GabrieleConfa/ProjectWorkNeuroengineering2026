function geom3 = addActiveSites(geom3, sites, params)
%ADDACTIVESITES Insert the 48 active sites (and the shank, if modelled) in the 3D geometry.
%   geom3 = addActiveSites(geom3, sites, params)
%
%   Owner: WP2. Status: SKELETON.
%
%   Inputs
%     geom3   struct from buildNerveGeometry3D
%     sites   table from loadActiveSites (centre_m, n, w, thickness_m, sideZ_m, sideXY_m)
%     params  struct from default_params() (params.electrode.siteModel, withShank)
%
%   Output: geom3 with, in addition,
%     siteFaces    (48,1) cell: face IDs through which each AS injects current
%     siteArea_m2  (48,1) area of those faces [m^2]
%     cellTissue   updated (code 5 = electrode, if sites are volumes)
%
%   Each AS is a box: multicuboid(sideXY, thickness, sideZ) rotated so that its
%   thickness is along n and its sideXY along w, then translated to centre_m.
%   Checked 2026-11-08: addCell and addVoid accept a 5 um box inside an extruded
%   cell. Not yet checked: cuff sites touching the nerve surface; AS crossing a
%   perineurium (r = 0.3 mm sites with patient P00).

arguments
    geom3 (1,1) struct
    sites table
    params (1,1) struct
end
error('project:notImplemented', 'addActiveSites: not implemented yet (WP2).');
end

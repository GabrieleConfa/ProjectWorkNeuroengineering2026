function geom3 = buildNerveGeometry3D(xsec, params)
%BUILDNERVEGEOMETRY3D Extrude the 2D cross-section into the 3D nerve + surrounding medium.
%   geom3 = buildNerveGeometry3D(xsec, params)
%
%   Owner: WP2. Status: SKELETON.
%
%   Inputs
%     xsec    struct from buildCrossSection2D
%     params  struct from default_params() (params.domain.length_m)
%
%   Output: struct
%     g             fegeometry object, 3D, spanning z = -L/2 ... +L/2 (device centred at z = 0)
%     cellTissue    (nCells,1) tissue code of each cell (1 endo ... 4 medium)
%     cellFascicle  (nCells,1) fascicle index, 0 outside
%
%   Checked 2026-11-08: extrude(fegeometry(dl), L) works; after extrude the
%   geometry starts at z = 0, so translate by -L/2. Cell numbering after extrude
%   must be mapped back to the 2D faces (do not assume the same order).

arguments
    xsec (1,1) struct
    params (1,1) struct
end
error('project:notImplemented', 'buildNerveGeometry3D: not implemented yet (WP2).');
end

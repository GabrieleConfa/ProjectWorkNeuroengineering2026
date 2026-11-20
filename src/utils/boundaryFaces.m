function tri = boundaryFaces(tets)
%BOUNDARYFACES Surface triangles of a set of linear tetrahedra (faces that belong to one element only).
%   tri = boundaryFaces(tets)
%
%   A face shared by two tetrahedra of the set is inside the region; a face that appears once is on
%   its boundary. Used to draw the surface of a region (a fascicle, a contact, the nerve) in 3D.
%   If the set is a slab of a longer region, the cut faces at the ends are boundary faces too, so
%   the drawing shows the region's cross-section there.
%
%   Input   tets  4 x nElem node indices (use the ORIGINAL node numbering: if the perineurium
%                 thin layer duplicated some nodes, map the copies back first, otherwise the faces
%                 between a fascicle and the epineurium do not match)
%   Output  tri   3 x nTri node indices of the boundary triangles (orientation not consistent)

arguments
    tets (4,:) double
end
F = [tets([1 2 3], :), tets([1 2 4], :), tets([1 3 4], :), tets([2 3 4], :)];
[~, ia, ic] = unique(sort(F, 1)', 'rows');
cnt = accumarray(ic, 1);
tri = F(:, ia(cnt == 1));
end

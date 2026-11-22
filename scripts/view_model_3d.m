%% VIEW_MODEL_3D  Interactive 3D view of a solved 3D model: geometry, mesh, potential.
%
% How to use
%   1. Run main_pipeline_3d first (it writes results/main3d/<id>_3d.mat).
%   2. Set the options below and press Run (or type: setup_paths; view_model_3d).
%   3. Rotate with the mouse (the rotate tool is switched on), zoom with the wheel, pan with the
%      hand tool of the figure toolbar.
%
% Windows
%   Figure 1 - geometry: nerve surface (beige, transparent), fascicles (one colour each, labelled
%              F1...F10 at the right end of the drawn range), contacts (solved ones red, the others
%              gold), cuff wall (grey, transparent). The cut ends of the drawn range show the
%              cross-section.
%   Figure 2 - potential of one solved contact, for 1 uA, on two planes: the cross-section through
%              the contact and the longitudinal plane through the nerve axis and the contact.
%              White = outside the model (the cuff wall is a void).
%   Both are also saved as PNG in results/main3d/.
%
% The nerve is 70 mm long and 3.4 mm wide: drawn whole, it is a thin stick. Draw a slab around the
% contacts you want to inspect (zRange_mm), or use zScale < 1 to compress z.

clear; close all;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);  setup_paths;

%% ===== OPTIONS ========================================================================
patientId     = "UA01";
zRange_mm     = [-13 -7.5];     % slab of nerve to draw [mm]; whole model: [-35 35]
zScale        = 1;              % 1 = true proportions; e.g. 0.2 compresses z five times
showNerve     = true;
showFascicles = true;
showCuff      = true;
contactsShown = "all";          % "all" | "solved" | "none"
showMeshEdges = false;          % draw the surface mesh (useful to inspect it; slow on long slabs)
potentialSite = 2;              % site id (one of the solved contacts) for figure 2; [] = no figure 2
potentialRes  = 241;            % grid points per side of the potential planes
planesOutside = false;          % also colour the muscle outside the nerve on the planes (it is
                                % ~0 V inside the cuff zone, where the insulating wall separates it)

%% ===== LOAD ===========================================================================
mm = 1e-3;   uA = 1e-6;
S = load(fullfile(root, 'results', 'main3d', patientId + "_3d.mat"), 'out');   out = S.out;
X  = double(out.mesh.nodes_m);   T = double(out.mesh.tets);
nO = out.mesh.nOriginalNodes;    dupOf = double(out.mesh.dupOf);
Torig = T;   isCopy = T > nO;   Torig(isCopy) = dupOf(T(isCopy) - nO);   % copies -> originals
elemCell   = double(out.mesh.elemCell);
elemTissue = out.cells.tissue(elemCell)';
elemFasc   = out.cells.fascicle(elemCell)';
elemSite   = out.cells.site(elemCell)';
patient = out.patient;   R = patient.nerveRadius_m;   F = patient.fascicles;   nF = height(F);
if isfield(patient, 'nerveSemiAxes_m'), nAx = patient.nerveSemiAxes_m;   rotN = patient.nerveRotation_deg;
else, nAx = [R R];   rotN = 0;                        % results saved before the oval (2026-11-22)
end
R = max(nAx);                                         % half-size of the drawn planes
p = out.params;   sites = out.sites;
zc = mean(reshape(X(3, T), 4, []), 1);                 % element centroids along z
inSlab = zc >= zRange_mm(1)*mm & zc <= zRange_mm(2)*mm;
fprintf('Drawing z = %g ... %g mm: %d of %d elements\n', zRange_mm, nnz(inSlab), numel(inSlab));

%% ===== FIGURE 1: GEOMETRY =============================================================
f1 = figure(Color="w", Position=[40 60 1200 800], Name="3D model - geometry");
ax = axes(f1);   hold(ax, "on");
edgeCol = "none";   if showMeshEdges, edgeCol = [0.2 0.2 0.2]; end
% Surface of a region = boundary faces of ALL its elements (whole length), then only the faces
% inside the slab are drawn: the tubes end cleanly at the slab limits instead of showing the
% jagged faces of the cut elements.
zLim = zRange_mm * mm;
inZ = @(tri) all(reshape(X(3, tri), 3, []) >= zLim(1) & reshape(X(3, tri), 3, []) <= zLim(2), 1);
slabFaces = @(sel) keepCols(boundaryFaces(Torig(:, sel)), inZ);
drawRegion = @(sel, col, alpha) patch(ax, Faces=slabFaces(sel)', Vertices=X'/mm, ...
    FaceColor=col, FaceAlpha=alpha, EdgeColor=edgeCol, EdgeAlpha=0.3, LineWidth=0.1, ...
    BackFaceLighting="lit", FaceLighting="gouraud");
if showNerve
    drawRegion(ismember(elemTissue, [1 3]), [0.98 0.85 0.62], 0.12);
end
cols = lines(nF);
if showFascicles
    for i = 1:nF
        sel = elemTissue == 1 & elemFasc == i;
        drawRegion(sel, cols(i,:), 0.35);
        th = linspace(0, 2*pi, 120);                    % outline of the fascicle at both slab ends
        for zEnd = zRange_mm
            plot3(ax, F.x_m(i)/mm + F.rEndo_m(i)/mm*cos(th), F.y_m(i)/mm + F.rEndo_m(i)/mm*sin(th), ...
                  zEnd*ones(size(th)), '-', Color=cols(i,:)*0.7, LineWidth=1);
        end
        text(ax, F.x_m(i)/mm, F.y_m(i)/mm, zRange_mm(2) + 0.3, "F" + i, FontWeight="bold", ...
             Color=cols(i,:)*0.7, HorizontalAlignment="center");
    end
end
switch contactsShown
    case "all",    showIds = sites.id;
    case "solved", showIds = out.solvedIds;
    otherwise,     showIds = [];
end
for s = showIds(:)'
    sel = elemSite == s & inSlab;
    if ~any(sel), continue; end
    col = [0.85 0.65 0.13];   if ismember(s, out.solvedIds), col = [0.85 0 0]; end
    drawRegion(sel, col, 1);
    c = sites.centre_m(sites.id == s, :) / mm;
    text(ax, c(1)*1.25, c(2)*1.25, c(3), "AS" + s, FontSize=7, Color=col*0.8);
end
if showCuff                                             % analytic tube: the wall is not meshed
    zc1 = max(zRange_mm(1), -p.cuff.length_m/2/mm);   zc2 = min(zRange_mm(2), p.cuff.length_m/2/mm);
    if zc2 > zc1
        co = nerveOutline(nAx + p.cuff.wallThickness_m, rotN, 90) / mm;   co = co([1:end 1], :);
        ni = nerveOutline(nAx, rotN, 90) / mm;                            ni = ni([1:end 1], :);
        surf(ax, [co(:,1) co(:,1)]', [co(:,2) co(:,2)]', repmat([zc1; zc2], 1, size(co, 1)), ...
             FaceColor=[0.3 0.3 0.3], FaceAlpha=0.08, EdgeColor="none");
        for zEnd = [zc1 zc2]                            % rings at the ends of the drawn wall
            plot3(ax, co(:,1), co(:,2), zEnd*ones(size(co, 1), 1), '-', Color=[0.3 0.3 0.3]);
            plot3(ax, ni(:,1), ni(:,2), zEnd*ones(size(ni, 1), 1), '-', Color=[0.3 0.3 0.3]);
        end
    end
end
finishAxes(ax, zScale);
title(ax, sprintf('Patient %s: z = %g ... %g mm (solved contacts red, others gold; cuff grey)', ...
      patientId, zRange_mm), FontSize=10);
exportgraphics(f1, fullfile(root, 'results', 'main3d', patientId + "_3d_view_geometry.png"), Resolution=150);

%% ===== FIGURE 2: POTENTIAL ON TWO PLANES ==============================================
if ~isempty(potentialSite)
    j = find(out.solvedIds == potentialSite);
    assert(~isempty(j), 'Site %d was not solved; solved: %s.', potentialSite, mat2str(out.solvedIds'));
    row = find(sites.id == potentialSite);
    tr = triangulation(T', X');
    u = double(out.V_per_A(:, j)) * uA;                 % volts for 1 uA
    c = sites.centre_m(row, :);
    er = [cosd(sites.theta_deg(row)), sind(sites.theta_deg(row)), 0];   % radial direction of the contact
    n = potentialRes;
    % cross-section through the contact
    [xg, yg] = meshgrid(linspace(-1.25, 1.25, n) * R);
    Vx = reshape(interpolateP1(tr, u, [xg(:), yg(:), c(3)*ones(n^2, 1)]), n, n);
    % longitudinal plane through the axis and the contact
    [sg, zg] = meshgrid(linspace(-1.25, 1.25, n) * R, linspace(zRange_mm(1), zRange_mm(2), 2*n) * mm);
    Vl = reshape(interpolateP1(tr, u, [sg(:)*er(1), sg(:)*er(2), zg(:)]), size(sg));
    if ~planesOutside                                   % nerve only
        Vx(reshape(ellipseQ([xg(:) yg(:)], nAx, rotN), size(xg)) > 1) = NaN;
        Vl(abs(sg) > 1 / sqrt(ellipseQ(er(1:2), nAx, rotN))) = NaN;     % nerve boundary along er
    end
    lim = prctile([Vx(:); Vl(:)], [1 99.5]) / 1e-3;

    f2 = figure(Color="w", Position=[80 40 1200 800], Name="3D model - potential");
    ax2 = axes(f2);   hold(ax2, "on");
    surf(ax2, xg/mm, yg/mm, c(3)/mm*ones(n), Vx/1e-3, EdgeColor="none", FaceColor="interp");
    surf(ax2, sg*er(1)/mm, sg*er(2)/mm, zg/mm, Vl/1e-3, EdgeColor="none", FaceColor="interp");
    for i = 1:nF                                        % fascicle outlines on the cross-section
        th = linspace(0, 2*pi, 120);
        plot3(ax2, (F.x_m(i) + F.rEndo_m(i)*cos(th))/mm, (F.y_m(i) + F.rEndo_m(i)*sin(th))/mm, ...
              c(3)/mm*ones(size(th)), 'k-', LineWidth=0.6);
    end
    sel = elemSite == potentialSite;
    patch(ax2, Faces=boundaryFaces(Torig(:, sel))', Vertices=X'/mm, FaceColor=[0.85 0 0], EdgeColor="none");
    if showFascicles                                    % faint fascicle surfaces for orientation
        for i = 1:nF
            s2 = elemTissue == 1 & elemFasc == i;
            patch(ax2, Faces=slabFaces(s2)', Vertices=X'/mm, FaceColor=cols(i,:), ...
                  FaceAlpha=0.07, EdgeColor="none");
        end
    end
    colormap(ax2, turbo);   clim(ax2, lim);   cb = colorbar(ax2);   cb.Label.String = 'V [mV] for 1 uA';
    finishAxes(ax2, zScale);
    title(ax2, sprintf('AS %d (%s, %s): potential for 1 uA on the cross-section and the longitudinal plane', ...
          potentialSite, string(sites.type(row)), out.results.location(j)), FontSize=10);
    exportgraphics(f2, fullfile(root, 'results', 'main3d', sprintf('%s_3d_view_AS%d.png', patientId, potentialSite)), ...
                   Resolution=150);
end

function A = keepCols(A, test)
%KEEPCOLS Columns of A for which test(A) is true (test returns a logical row).
A = A(:, test(A));
end

function finishAxes(ax, zScale)
%FINISHAXES Proportions, labels, light, 3D view and the rotate tool.
daspect(ax, [1 1 1/zScale]);   axis(ax, "tight");   grid(ax, "on");   box(ax, "on");
xlabel(ax, 'x [mm]');   ylabel(ax, 'y [mm]');   zlabel(ax, 'z [mm] (nerve axis)');
view(ax, [-38 22]);   camlight(ax, "headlight");   lighting(ax, "gouraud");
rotate3d(ax, "on");
end

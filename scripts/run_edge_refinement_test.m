function T = run_edge_refinement_test(edgeSizes_um, padIds) % Does refining the pad edges stop the drift of the access resistance?
%% 0. Idea being tested
% Current crowds on the edges and corners of a flat pad. If the mesh is too coarse there, the crowding is smoothed,
% the pad behaves as if it were wider, and the access resistance R is underestimated. Here the production mesh is kept
% (size factor 1) and ONLY the element size on the pad edges (p.mesh3d.hEdge_m) is reduced. If R then converges, the
% edge explanation is supported; if not, another cause is at work. Nothing in results/main3d is overwritten.
if nargin<1, edgeSizes_um=[0 20 10]; end % 0 = no edge refinement (production mesh).
if nargin<2, padIds=[1 5 17]; end % AS1 internal pad, AS5 and AS17 cuff pads.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root); setup_paths(); % Resolve all sources relative to this file.
outDir=fullfile(root,'results','convergence'); if ~isfolder(outDir),mkdir(outDir);end % Output folder.
%% 1. Same personal case as start_nerve
previous=getenv('NERVE_CASE'); restore=onCleanup(@()setenv('NERVE_CASE',previous)); setenv('NERVE_CASE','personal'); %#ok<NASGU> % Select the personal geometry.
[p,patient,sites]=personalDevice3D(pw6_3d_params(),quicklook3DParams()); % Same anatomy, pads and cuff.
p.solve3d.solver="pcg"; % Same solver as the production run.
%% 2. Mesh and solve for each edge size
n=numel(edgeSizes_um); R=nan(n,numel(padIds)); N=nan(n,1); sec=N; % Result arrays.
for k=1:n % One mesh per edge size.
    q=p; q.mesh3d.hEdge_m=edgeSizes_um(k)*1e-6; q.mesh3d.hEdgePads=padIds; % Edge element size in metres (0 = off), only for the pads under test to keep the mesh small.
    fprintf('Edge size %g um\n',edgeSizes_um(k)); % Progress report.
    r=solvePersonalPads(q,patient,sites,padIds); % Mesh, assemble and solve the chosen pads.
    R(k,:)=r.R_kOhm; N(k)=r.nTets; sec(k)=r.time_s; % Access resistances and mesh size.
    fprintf('   %d tetrahedra, R = %s kOhm\n',N(k),mat2str(round(r.R_kOhm,3))); % Report.
end
%% 3. Table and figure
T=table(edgeSizes_um(:),N,R(:,1),R(:,2),R(:,3),'VariableNames',{'edgeSize_um','tetrahedra',sprintf('R_AS%d_kOhm',padIds(1)),sprintf('R_AS%d_kOhm',padIds(2)),sprintf('R_AS%d_kOhm',padIds(3))}); % Summary table.
writetable(T,fullfile(outDir,'edge_refinement.csv')); disp(T); % Save and print.
fig=figure(Visible="off",Color="w",Position=[60 60 1000 520],Theme="light"); ax=axes(fig); hold(ax,'on'); grid(ax,'on'); % Plot R against the edge size.
x=1:n; for j=1:numel(padIds), plot(ax,x,R(:,j)/R(1,j)*100,'-o',LineWidth=1.6,DisplayName=sprintf('AS%d',padIds(j))); end % R relative to the production mesh.
set(ax,XTick=x,XTickLabel=compose('%g µm\n%d tets',edgeSizes_um(:),N)); % Label the meshes.
ylabel(ax,'Access resistance (% of production mesh)'); xlabel(ax,'Element size on pad edges'); legend(ax,Location='best'); % Labels.
title(ax,'Effect of refining only the pad edges (0 µm = production mesh)'); exportgraphics(fig,fullfile(outDir,'edge_refinement.png'),Resolution=170); close(fig); % Save figure.
end

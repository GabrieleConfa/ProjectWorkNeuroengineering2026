function summary = run_convergence_study(scales, padIds) % Mesh-convergence study of the personal 60 mm model.
%% 0. What this study answers
% The electrical problem: a current I = 1 uA leaves one active pad and flows through anisotropic
% endoneurium, a thin perineurium layer (potential jump), epineurium and muscle, to the grounded
% outer boundary. The potential V solves div(sigma grad V) = 0 (quasi-static volume conductor).
% The finite-element solution depends on the tetrahedra size h. Here the SAME geometry is meshed
% with 3 (or more) element sizes; if the quantities of interest stop changing as h shrinks, the
% mesh is fine enough for them. Every length in the mesh controls is multiplied by one factor.
% Nothing in the anatomy, leads or cuff is changed, and results/main3d is not overwritten.
if nargin<1, scales=[1.96 1.4 1.0 0.714]; end % Size factors relative to the production mesh (ratio 1.4 between levels; 0.714 is finer than production).
if nargin<2, padIds=[1 5 17]; end % AS1 = internal pad, AS5 and AS17 = cuff pads (AS17 has the largest indentation).
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root); setup_paths(); % Resolve all sources relative to this file.
outDir=fullfile(root,'results','convergence'); if ~isfolder(outDir),mkdir(outDir);end % Folder for figures and tables.
%% 1. Build the same personal case used by start_nerve
previous=getenv('NERVE_CASE'); restore=onCleanup(@()setenv('NERVE_CASE',previous)); setenv('NERVE_CASE','personal'); %#ok<NASGU> % Select the personal geometry.
[p,patient,sites]=personalDevice3D(pw6_3d_params(),quicklook3DParams()); % Same anatomy, pads and cuff as the production run.
p.solve3d.solver="pcg"; % Memory-efficient iterative solver, as in the production run.
base=p.mesh3d; names=["hContact_m","hEndo_m","hNerve_m","hCoarse_m","hMax_m"]; % Size controls that are scaled together.
%% 2. Mesh and solve at every refinement level
nL=numel(scales); nP=numel(padIds); % Number of meshes and of pads.
R=nan(nL,nP); spread=R; V0=R; Ezmax=R; AF=R; cons=R; N=nan(nL,1); nNodes=N; tMesh=N; % Result arrays.
profile=cell(nL,nP); zAxis=cell(nL,nP); % Axial potential along the probe line.
for k=1:nL % Loop from the coarsest to the finest mesh.
    q=p; for n=names, q.mesh3d.(n)=base.(n)*scales(k); end % Scale all element-size controls by the same factor.
    fprintf('Level %d/%d: size factor %.2f\n',k,nL,scales(k)); % Progress report.
    r=solvePersonalPads(q,patient,sites,padIds); % Mesh, assemble and solve the selected pads.
    N(k)=r.nTets; nNodes(k)=r.nNodes; tMesh(k)=r.time_s; % Mesh statistics.
    R(k,:)=r.R_kOhm; spread(k,:)=r.spread_pct; V0(k,:)=r.Vprobe_V; Ezmax(k,:)=r.EzMax_Vpm; AF(k,:)=r.AF_Vpm2; cons(k,:)=r.cons_pct; % Store the metrics.
    profile(k,:)=r.profile; zAxis(k,:)=r.z; % Store the profiles.
    fprintf('   %d tetrahedra, R = %s kOhm\n',N(k),mat2str(round(r.R_kOhm,3))); % Report access resistances.
end
%% 3. Differences, observed order and Richardson extrapolation
metrics=struct('R_access_kOhm',R,'V_probe_V',V0,'peak_Ez_Vpm',Ezmax,'activating_fn_Vpm2',AF); % Quantities compared.
fn=fieldnames(metrics); rows=table(); r=scales(end-2)/scales(end-1); % Refinement ratio between consecutive levels (order is estimated on the three finest).
for i=1:numel(fn) % Loop over quantities.
    for j=1:nP % Loop over pads.
        f=metrics.(fn{i})(:,j); % Values on the meshes, coarse to fine.
        change=100*abs(diff(f))./abs(f(2:end)); % Relative change between consecutive meshes (%).
        f3=f(end-2:end); % Three finest meshes.
        order=NaN; ext=NaN; gci=NaN; % Defaults for non-monotone sequences.
        if nL>=3 && (f3(1)-f3(2))*(f3(2)-f3(3))>0 % Monotone convergence is required for an order estimate.
            order=log((f3(1)-f3(2))/(f3(2)-f3(3)))/log(r); % Observed order of convergence p.
            ext=f3(3)+(f3(3)-f3(2))/(r^order-1); % Richardson extrapolation to h -> 0.
            gci=100*1.25*abs((f3(3)-f3(2))/f3(3))/(r^order-1); % Grid convergence index (safety factor 1.25).
            if order<0.5 || order>2.5, ext=NaN; gci=NaN; end % Linear tetrahedra cannot converge faster than ~2nd order: outside this range the extrapolation is not trustworthy.
        end
        rows=[rows;table(string(fn{i}),padIds(j),f(end),change(end),order,ext,gci,'VariableNames',{'quantity','pad','finest','changeLastStep_pct','observedOrder','extrapolated','GCI_pct'})]; %#ok<AGROW> % One row per quantity and pad.
    end
end
summary=rows; writetable(rows,fullfile(outDir,'convergence_summary.csv')); % Numbers for the report.
levels=table((1:nL)',scales(:),N,nNodes,tMesh,'VariableNames',{'level','sizeFactor','tetrahedra','nodes','seconds'}); writetable(levels,fullfile(outDir,'convergence_levels.csv')); % Mesh sizes.
% Consistency check against the production solution, which was computed on the same finest mesh.
prod=fullfile(root,'results','main3d','Cuff60_3d.mat'); check=nan(1,nP); % Saved result of start_nerve.
if isfile(prod) && any(abs(scales-1)<1e-12)
    S=load(prod,'out'); [~,ix]=ismember(padIds,S.out.results.id); check=S.out.results.R_kOhm(ix)'; % Saved access resistances.
end
save(fullfile(outDir,'convergence.mat'),'scales','padIds','N','R','spread','V0','Ezmax','AF','cons','profile','zAxis','rows','levels','check'); % Keep everything.
%% 4. Figures
plotConvergence(outDir,scales,padIds,N,metrics,profile,zAxis,rows,check); % Summary and profile figures.
disp(levels); disp(rows); % Print the tables.
end
%% Local functions
function plotConvergence(outDir,scales,padIds,N,metrics,profile,zAxis,rows,check) % Explanatory figures.
nL=numel(scales); nP=numel(padIds); cols=lines(nL); lab=compose('%d tets (h x %.2f)',N,scales(:)); % Colours and legend labels.
kind=@(id)ternary(id==1,'internal pad','cuff pad'); % Pad type for titles.
%% Figure 1: values versus mesh size, with extrapolation
fig=figure(Visible="off",Color="w",Position=[40 40 1250 800],Theme="light"); tl=tiledlayout(fig,2,2,TileSpacing="compact",Padding="compact"); % Four quantities.
labels={'Access resistance (kΩ)','Potential on probe line at pad plane (V per A)','Peak |E_z| on probe line (V/m per A)','Activating function d²V/dz² (V/m² per A)'}; % Axis labels.
fn=fieldnames(metrics); mk='osd^'; % One marker per pad.
for i=1:4 % One panel per quantity.
    ax=nexttile(tl); hold(ax,'on'); grid(ax,'on'); % Axes.
    for j=1:nP % One line per pad, normalized to its own finest value.
        f=metrics.(fn{i})(:,j); plot(ax,N,100*f/f(end),['-' mk(j)],LineWidth=1.6,DisplayName=sprintf('AS%d (%s)',padIds(j),kind(padIds(j)))); % Percentage of the finest-mesh value.
    end
    set(ax,XScale='log'); yline(ax,100,'k:','HandleVisibility','off'); % Reference: finest mesh = 100 %.
    if i==3, ylim(ax,[50 200]); text(ax,0.02,0.95,'values above 200 % clipped (coarse-mesh spike)','Units','normalized','VerticalAlignment','top'); end % Peak field is spiky on coarse P1 meshes.
    xlabel(ax,'Number of tetrahedra'); ylabel(ax,'Value as % of finest mesh'); title(ax,labels{i}); % Labels.
    if i==1, legend(ax,Location='best'); end % One legend.
end
title(tl,'Mesh convergence: each curve should flatten near 100 % as the mesh is refined (production mesh = 3rd point)'); exportgraphics(fig,fullfile(outDir,'convergence_values.png'),Resolution=170); close(fig); % Save.
%% Figure 2: axial profiles on the three meshes
fig=figure(Visible="off",Color="w",Position=[40 40 1250 850],Theme="light"); tl=tiledlayout(fig,nP,2,TileSpacing="compact",Padding="compact"); % Pads x (V, Ez).
for j=1:nP
    ax=nexttile(tl); hold(ax,'on'); grid(ax,'on'); % Potential panel.
    for k=1:nL, plot(ax,zAxis{k,j},profile{k,j},Color=cols(k,:),LineWidth=1.5,DisplayName=lab{k}); end % One curve per mesh.
    xlabel(ax,'z from pad centre (mm)'); ylabel(ax,'V per A (V/A)'); title(ax,sprintf('AS%d (%s): potential on probe line',padIds(j),kind(padIds(j)))); if j==1,legend(ax,Location='best'); end % Labels.
    ax=nexttile(tl); hold(ax,'on'); grid(ax,'on'); % Axial field panel.
    for k=1:nL, plot(ax,zAxis{k,j},-gradient(profile{k,j},zAxis{k,j}*1e-3),Color=cols(k,:),LineWidth=1.5); end % E_z = -dV/dz.
    xlabel(ax,'z from pad centre (mm)'); ylabel(ax,'E_z per A (V/m per A)'); title(ax,sprintf('AS%d: axial electric field',padIds(j))); % Labels.
end
title(tl,'Axial profiles 0.3 mm from the active face, on each mesh'); exportgraphics(fig,fullfile(outDir,'convergence_profiles.png'),Resolution=170); close(fig); % Save.
%% Figure 3: change in the last refinement step
fig=figure(Visible="off",Color="w",Position=[40 40 1000 520],Theme="light"); ax=axes(fig); % Bar chart.
names=["R access","V probe","peak Ez","Act. fn"]; C=zeros(numel(fn),nP); % Matrix quantities x pads.
for i=1:numel(fn), for j=1:nP, C(i,j)=rows.changeLastStep_pct(rows.quantity==string(fn{i}) & rows.pad==padIds(j)); end, end % Last-step changes.
b=bar(ax,C); for j=1:nP, b(j).DisplayName=sprintf('AS%d (%s)',padIds(j),kind(padIds(j))); end % One bar group per quantity.
set(ax,XTickLabel=names); ylabel(ax,'Change between the two finest meshes (%)'); grid(ax,'on'); legend(ax,Location='northwest'); % Labels.
title(ax,'How much each quantity still changes in the last refinement step'); exportgraphics(fig,fullfile(outDir,'convergence_change.png'),Resolution=170); close(fig); % Save.
if all(isfinite(check)), fprintf('Production check (saved results vs this study, scale 1 mesh): max difference %.3g %%\n',max(100*abs(metrics.R_access_kOhm(abs(scales-1)<1e-12,:)-check)./check)); end % Compare with start_nerve output.
end
function v=ternary(c,a,b) % Small inline conditional.
if c, v=a; else, v=b; end
end

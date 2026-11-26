function T = run_pdms_sensitivity(padIds) % How much do the PDMS conductivity and the sleeve length matter? Where does the current go?
%% Purpose
% (1) PDMS conductivity: unfilled PDMS is an insulator, but filled (conductive) PDMS exists. The lead bodies and the sleeve are
%     meshed once per case and solved with different conductivities.
% (2) Sleeve length: the oval PDMS sleeve is assumed 50 mm long (cuff length of the group). A shorter sleeve leaves more nerve
%     surface in contact with the muscle.
% (3) Current paths: share of the injected current that reaches the grounded end planes of the nerve (z = +/-30 mm) and the
%     grounded outer muscle surface.
% Nothing in results/main3d is overwritten.
if nargin<1, padIds=[1 5 17]; end % Internal pad AS1, cuff pads AS5 and AS17.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root); setup_paths(); outDir=fullfile(root,'results','convergence'); if ~isfolder(outDir),mkdir(outDir);end
previous=getenv('NERVE_CASE'); restore=onCleanup(@()setenv('NERVE_CASE',previous)); setenv('NERVE_CASE','personal'); %#ok<NASGU>
cases={'PDMS insulator (baseline), sleeve 50 mm',1e-7,50e-3;'PDMS conductive 1e-3 S/m, sleeve 50 mm',1e-3,50e-3;'PDMS conductive 1e-1 S/m, sleeve 50 mm',1e-1,50e-3; ...
       'PDMS insulator, sleeve 36 mm',1e-7,36e-3;'PDMS insulator, sleeve 22 mm',1e-7,22e-3}; % Sleeves must cover the pads (|z| <= 10.4 mm).
rows=cell(size(cases,1),1);
for k=1:size(cases,1)
    [p,patient,sites]=personalDevice3D(pw6_3d_params(),quicklook3DParams()); p.solve3d.solver="pcg"; % Same case as the production run.
    p.electrode.pdmsSigma_Spm=cases{k,2}; p.cuff.length_m=cases{k,3}; % Case parameters.
    fprintf('%s\n',cases{k,1}); r=solvePersonalPads(p,patient,sites,padIds); % Mesh and solve the chosen pads.
    rows{k}=table(string(cases{k,1}),cases{k,2},1e3*cases{k,3},r.R_kOhm(1),r.R_kOhm(2),r.R_kOhm(3),mean(r.endPlanes_pct),mean(r.sideMuscle_pct),'VariableNames', ...
        {'case','pdmsSigma_Spm','sleeve_mm',sprintf('R_AS%d_kOhm',padIds(1)),sprintf('R_AS%d_kOhm',padIds(2)),sprintf('R_AS%d_kOhm',padIds(3)),'current_endPlanes_pct','current_muscleSide_pct'});
end
T=vertcat(rows{:}); writetable(T,fullfile(outDir,'pdms_sensitivity.csv')); disp(T);
fig=figure(Visible="off",Color="w",Position=[60 60 1100 480],Theme="light"); tl=tiledlayout(fig,1,2,TileSpacing="compact"); % Figure.
ax=nexttile(tl); bar(ax,categorical(T.case,T.case),T{:,4:6}); ylabel(ax,'Access resistance (k\Omega)'); legend(ax,compose('AS%d',padIds),Location='northeast'); grid(ax,'on'); title(ax,'Access resistance'); ax.XTickLabelRotation=20; % Resistances.
ax=nexttile(tl); bar(ax,categorical(T.case,T.case),[T.current_endPlanes_pct T.current_muscleSide_pct],'stacked'); ylabel(ax,'Share of the injected current (%)'); legend(ax,{'nerve end planes','muscle surface'},Location='best'); grid(ax,'on'); title(ax,'Where the current goes'); ax.XTickLabelRotation=20; % Current paths.
exportgraphics(fig,fullfile(outDir,'pdms_sensitivity.png'),Resolution=170); close(fig);
end

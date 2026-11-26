%% RUN_CONVERGENCE_3D  Mesh-convergence study of the main 3D model (main_pipeline_3d).
%
% What it answers
%   The FE solution depends on the element size h. The same geometry is meshed with all element-size
%   targets multiplied by one factor ("level"), the same contacts are solved on every mesh, and the
%   same quantities are evaluated AT THE SAME PHYSICAL POINTS. If they stop changing as h shrinks,
%   the production mesh (factor 1) is fine enough for them, and the change gives its error.
%
% Method (docs/design.md §4, proposal of 2026-11-23, and the algorithm of Marco's
% scripts/run_convergence_study.m on branch cuff-oval-viewer-convergence, merged 2026-11-26)
%   1. Levels: size factors with a constant ratio r = 1.3 between consecutive meshes, from coarse to
%      fine (the finest is finer than production). The effective size of an unstructured mesh is
%      h = (volume / number of elements)^(1/3); the ratio of the factors is used as r.
%   2. Remeshing noise: our meshes are not nested; every level is a new tetrahedralisation. The
%      production level is meshed a second time with all sizes x 1.02 (a 2 % change that should not
%      change the answer): the difference is the noise floor. A change between levels below it says
%      nothing.
%   3. Quantities, for each test contact: access resistance R; V at the centre of the nearest
%      fascicle in the contact plane; V at the centre of the farthest fascicle; activating function
%      (second difference over one internode, p.plot3d.afSpacing_m) on the axis of the nearest
%      fascicle at the contact plane. Per-mesh sanity check: the current through the ground equals the
%      injected current (conservErr_pct).
%   4. Error estimate, Richardson in the GCI form with checks (Roache; Celik et al. 2008, J Fluids
%      Eng 130:078001, not yet in Papers/): from the three finest levels, the observed order
%      p = ln((f1 - f2)/(f2 - f3)) / ln r. If the changes are monotonic and 0.5 <= p <= 3, report the
%      extrapolated value and GCI = 1.25 |(f3 - f2)/f3| / (r^p - 1). Otherwise do not extrapolate:
%      report the change between the two finest levels (criterion of Davis 2023: < 0.5 % on the
%      threshold) and compare it with the noise floor.
%   5. Edge refinement (Marco's result on his Gmsh model: R rises by 2-5 % with 5 um elements on the
%      pad edges, then stops changing): one extra mesh at factor 1 with the edge refinement toggled
%      (switched OFF if production has it on, as since D-14; ON at edgeTest_m otherwise).
%   The mesher is the production one (p.mesh3d.mesher). Every level has its own size targets, hence
%   its own mesh key: the meshes are saved in results/mesh/ like the production one and never
%   overwrite it.
%
% Outputs (results/convergence3d/): convergence_levels.csv, convergence_summary.csv,
%   convergence_values.png, convergence_change.png, convergence.mat
% Run: setup_paths; run_convergence_3d
% Time / memory: the geometry is built once (minutes); then one mesh per level (minutes) and one pcg
%   solve per test contact. The finest level has ~2 M nodes: pcg for every level (the direct solver
%   would not fit in 32 GB there). Total ~1-1.5 h. Do not run it together with another big job.

clear;  close all;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);  setup_paths;
outDir = fullfile(root, 'results', 'convergence3d');
if ~exist(outDir, 'dir'), mkdir(outDir); end
mm = 1e-3;

%% Settings of the study ---------------------------------------------------------------------------
scales     = [1.69 1.3 1.0 0.77];   % size factors, coarse -> fine (ratio 1.3); 1.0 = production
noiseScale = 1.02;                  % second mesh of the production level, for the noise floor
edgeTest_m = 5e-6;                  % edge size of the toggled run when production has none (0 = skip the run)
testIds    = [1 2 5];               % AS1 (intraneural, F2), AS2 (intraneural, F7), AS5 (cuff pad)

p = pw6_3d_params();
p.solve3d.solver = "pcg";           % same solver on every level (memory)
p.solve3d.saveSingle = false;

%% Meshes and solves ---------------------------------------------------------------------------------
runs = [num2cell(scales), {noiseScale}];
kinds = [repmat("level", 1, numel(scales)), "noise"];
edgeAlt = edgeTest_m;   if p.mesh3d.hEdge_m > 0, edgeAlt = 0; end   % toggled edge size
if edgeTest_m > 0 || p.mesh3d.hEdge_m > 0, runs{end+1} = 1.0; kinds(end+1) = "edge"; end
nR = numel(runs);   nP = numel(testIds);
Q = struct('R_kOhm', nan(nR, nP), 'V_near_VpA', nan(nR, nP), 'V_far_VpA', nan(nR, nP), ...
           'AF_near_Vpm2pA', nan(nR, nP));
info = table((1:nR)', [runs{:}]', kinds', nan(nR, 1), nan(nR, 1), nan(nR, 1), nan(nR, 1), nan(nR, 1), ...
             'VariableNames', {'run', 'sizeFactor', 'kind', 'nodes', 'tetrahedra', 'h_eff_um', 'maxConservErr_pct', 'seconds'});
geom = [];
for k = 1:nR
    t0 = tic;
    q = p;
    if kinds(k) == "edge", q.mesh3d.hEdge_m = edgeAlt; end
    fprintf('\n=== Run %d/%d: %s, size factor %.2f%s ===\n', k, nR, kinds(k), runs{k}, ...
            string(repmat(sprintf(', edge refinement %.0f um (0 = off)', edgeAlt * 1e6), 1, kinds(k) == "edge")));
    M = build3DModel(q, struct('savePatient', false, 'sizeScale', runs{k}, 'geom', geom, 'verbose', k == 1));
    geom = M.geom;                                            % the geometry does not depend on h
    [V, res] = solve3DContacts(M, testIds(:), struct('verbose', true));
    tr = triangulation(M.tl.tets', M.tl.nodes_m');
    F = M.patient.fascicles;   S = M.sites;   dzAF = q.plot3d.afSpacing_m;
    for j = 1:nP
        s = S(S.id == testIds(j), :);
        d = hypot(F.x_m - s.centre_m(1), F.y_m - s.centre_m(2));
        [~, iN] = min(d);   [~, iF] = max(d);                 % nearest and farthest fascicle centres
        pts = [F.x_m(iN) * [1 1 1]', F.y_m(iN) * [1 1 1]', s.z_m + [-dzAF 0 dzAF]';
               F.x_m(iF), F.y_m(iF), s.z_m];
        v = interpolateP1(tr, double(V(:, j)), pts);
        Q.R_kOhm(k, j) = res.R_kOhm(j);
        Q.V_near_VpA(k, j) = v(2);   Q.V_far_VpA(k, j) = v(4);
        Q.AF_near_Vpm2pA(k, j) = (v(1) - 2 * v(2) + v(3)) / dzAF^2;
    end
    info.nodes(k) = M.nN;   info.tetrahedra(k) = M.nE;
    info.h_eff_um(k) = (sum(M.vol) / M.nE)^(1/3) * 1e6;
    info.maxConservErr_pct(k) = max(res.conservErr_pct);
    info.seconds(k) = toc(t0);
    fprintf('   %d nodes, %d tetrahedra, h_eff = %.0f um, R = %s kOhm  [%.0f s]\n', M.nN, M.nE, ...
            info.h_eff_um(k), mat2str(round(res.R_kOhm', 3)), info.seconds(k));
    clear M V tr
end
writetable(info, fullfile(outDir, 'convergence_levels.csv'));

%% Analysis: changes, noise floor, observed order, Richardson / GCI --------------------------------
iLev = find(kinds == "level");   iProd = iLev(abs([runs{iLev}] - 1) < 1e-12);
iNoise = find(kinds == "noise");   iEdge = find(kinds == "edge");
r = scales(end-2) / scales(end-1);                            % refinement ratio of the three finest
names = fieldnames(Q);   rows = table();
for i = 1:numel(names)
    for j = 1:nP
        f = Q.(names{i})(iLev, j);                             % coarse -> fine
        change = 100 * abs(diff(f)) ./ abs(f(2:end));
        noise = 100 * abs(Q.(names{i})(iNoise, j) - Q.(names{i})(iProd, j)) / abs(Q.(names{i})(iProd, j));
        edge = NaN;
        if ~isempty(iEdge)
            edge = 100 * (Q.(names{i})(iEdge, j) - Q.(names{i})(iProd, j)) / abs(Q.(names{i})(iProd, j));
        end
        f3 = f(end-2:end);   order = NaN;   ext = NaN;   gci = NaN;
        if (f3(1) - f3(2)) * (f3(2) - f3(3)) > 0               % monotonic on the three finest levels
            order = log((f3(1) - f3(2)) / (f3(2) - f3(3))) / log(r);
            if order >= 0.5 && order <= 3
                ext = f3(3) + (f3(3) - f3(2)) / (r^order - 1);
                gci = 100 * 1.25 * abs((f3(3) - f3(2)) / f3(3)) / (r^order - 1);
                verdict = "GCI valid";
            else
                verdict = "order outside 0.5-3: no extrapolation";
            end
        else
            verdict = "not monotonic: no extrapolation";
        end
        if change(end) <= noise, verdict = verdict + "; last change below the noise floor"; end
        rows = [rows; table(string(names{i}), testIds(j), f(end), f(iLev == iProd), change(end), noise, ...
                            order, ext, gci, edge, verdict, 'VariableNames', ...
                {'quantity', 'contact', 'finest', 'production', 'changeLastStep_pct', 'noiseFloor_pct', ...
                 'observedOrder', 'extrapolated', 'GCI_pct', 'edgeToggleEffect_pct', 'verdict'})]; %#ok<AGROW>
    end
end
writetable(rows, fullfile(outDir, 'convergence_summary.csv'));
save(fullfile(outDir, 'convergence.mat'), 'scales', 'runs', 'kinds', 'testIds', 'Q', 'info', 'rows');
disp(info);  disp(rows);

%% Figures ---------------------------------------------------------------------------------------------
% 1) every quantity vs number of tetrahedra, in % of its value on the finest level
labels = {'Access resistance R', 'V at the nearest fascicle centre', 'V at the farthest fascicle centre', ...
          'Activating function on the nearest fascicle axis'};
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [40 40 1250 800]);
tlo = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact');
mk = 'osd^';
for i = 1:numel(names)
    ax = nexttile(tlo);  hold(ax, 'on');  grid(ax, 'on');
    for j = 1:nP
        f = Q.(names{i})(iLev, j);
        plot(ax, info.tetrahedra(iLev), 100 * f / f(end), ['-' mk(j)], 'LineWidth', 1.5, ...
             'DisplayName', sprintf('AS%d', testIds(j)));
        plot(ax, info.tetrahedra(iNoise), 100 * Q.(names{i})(iNoise, j) / f(end), mk(j), ...
             'MarkerFaceColor', 'k', 'HandleVisibility', 'off');
        if ~isempty(iEdge)
            plot(ax, info.tetrahedra(iEdge), 100 * Q.(names{i})(iEdge, j) / f(end), mk(j), ...
                 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
        end
    end
    set(ax, 'XScale', 'log');  yline(ax, 100, 'k:', 'HandleVisibility', 'off');
    xlabel(ax, 'number of tetrahedra');  ylabel(ax, '% of the finest-level value');  title(ax, labels{i});
    if i == 1, legend(ax, 'Location', 'best'); end
end
title(tlo, {'Mesh convergence of the main 3D model: each curve should flatten near 100 % as the mesh is refined', ...
            sprintf('production = factor 1.0; black = remeshed at x1.02 (noise floor); red = production with edge size %g um (0 = off)', edgeAlt * 1e6)});
exportgraphics(fig, fullfile(outDir, 'convergence_values.png'), 'Resolution', 150);  close(fig);

% 2) change in the last refinement step vs the noise floor
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [40 40 1000 520]);
ax = axes(fig);
C = reshape(rows.changeLastStep_pct, nP, [])';   N = reshape(rows.noiseFloor_pct, nP, [])';
b = bar(ax, C);  hold(ax, 'on');
for j = 1:nP
    b(j).DisplayName = sprintf('AS%d', testIds(j));
    plot(ax, b(j).XEndPoints, N(:, j), 'k_', 'MarkerSize', 14, 'LineWidth', 2, 'HandleVisibility', 'off');
end
set(ax, 'XTickLabel', {'R', 'V near', 'V far', 'AF near'});
ylabel(ax, 'change between the two finest levels [%]');  grid(ax, 'on');  legend(ax, 'Location', 'northwest');
title(ax, 'How much each quantity still changes in the last refinement (black dash = remeshing noise floor)');
exportgraphics(fig, fullfile(outDir, 'convergence_change.png'), 'Resolution', 150);  close(fig);
fprintf('Wrote %s\n', outDir);

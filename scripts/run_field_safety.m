%% RUN_FIELD_SAFETY  Electric-field exposure of the tissue for the stimulation and block patterns.
%
% Why
%   The block study pushes the anodal current up to milliamperes. Near a 0.2 x 0.2 mm contact this
%   gives fields of hundreds to thousands of V/cm, the range where cell membranes electroporate:
%   the tissue is damaged, and its conductivity rises (a field-dependent sigma(|E|) that our linear
%   FEM ignores). Instead of building a nonlinear FEM (no superposition, one 3D Newton solve per
%   pattern and current), this script checks WHERE and AT WHICH CURRENT the field of the linear
%   model exceeds given levels. If a setting stays below the electroporation threshold, the linear
%   model is self-consistent there; if not, the setting is unsafe anyway.
%   Plan and reasons: docs/design.md §5.5 (roadmap item 5).
%
% What it does
%   For each pattern (the stimulation cathode and the anode configurations of
%   notebooks/block_anode_configurations.ipynb): fieldExposure gives the |E| distribution per tissue
%   for 1 A; then, for each field level of p.safety.fieldLevels_Vpm and a range of currents, the
%   volume of tissue above the level (exposedVolume), and the "critical current" at which that volume
%   reaches p.safety.exposedVolume_m3 (1e-3 mm^3 = a (100 um)^3 cube, [assumed]). The perineurium
%   (thin layer) is reported separately: the current at which the field across it exceeds each level
%   at any node.
%
% Outputs (results/safety/)
%   field_exposure_<id>.csv   critical currents per pattern, tissue and field level
%   field_exposure_<id>.png   exposed volume vs current, one panel per pattern
%
% Run: setup_paths; run_field_safety           (needs results/main3d/<id>_3d.mat)
% Time: ~1 min per pattern (field on 4.5 M elements). Memory: a few GB; do not run during a solve.

clear;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);  setup_paths;
mm = 1e-3;

p = pw6_3d_params();                       % safety levels (config/default_params.m, p.safety)
b = block_params();                        % patient id and the stimulation current context
S = load(fullfile(root, 'results', 'main3d', b.patientId + "_3d.mat"), 'out');   out = S.out;
outDir = fullfile(root, 'results', 'safety');
if ~exist(outDir, 'dir'), mkdir(outDir); end
levels = p.safety.fieldLevels_Vpm;
Vcrit = p.safety.exposedVolume_m3;
currents = logspace(log10(1e-6), log10(10e-3), 60)';         % 1 uA ... 10 mA (total pattern current)
fprintf('Field exposure of %s; solved contacts %s\n', b.patientId, mat2str(out.solvedIds(:)'));
if isnan(p.safety.electroporation_Vpm)
    fprintf('Electroporation threshold not set yet (p.safety.electroporation_Vpm = NaN): levels %s V/cm reported.\n', ...
            mat2str(levels / 100));
end

%% Patterns -----------------------------------------------------------------------------------------
% name, contact ids, weights (fraction of the total current; equal split for the anodes). The same
% list as ANODE_CONFIGS in notebooks/block_anode_configurations.ipynb, plus the stimulation cathode.
sites = out.sites;
F = out.patient.fascicles;   target = 7;                     % "nearest" = closest shank contact to F7
nearest = zeros(1, 4);
for sec = 1:4
    sh = sites(sites.section == sec & sites.type == "intraneural", :);
    [~, k] = min(hypot(sh.r_m .* cosd(sh.theta_deg) - F.x_m(target), sh.r_m .* sind(sh.theta_deg) - F.y_m(target)));
    nearest(sec) = sh.id(k);
end
eq = @(ids) ones(1, numel(ids)) / numel(ids);
pat = struct('name', {}, 'ids', {}, 'w', {});
pat(end+1) = struct('name', "stimulation cathode AS2", 'ids', 2, 'w', -1);
pat(end+1) = struct('name', "same contact AS44", 'ids', 44, 'w', 1);
pat(end+1) = struct('name', "station +4, all 6", 'ids', 43:48, 'w', eq(43:48));
pat(end+1) = struct('name', "station +4, 2 cuff pads", 'ids', [47 48], 'w', eq([47 48]));
cuff8 = [29 30 35 36 41 42 47 48];
shank16 = setdiff(25:48, cuff8);
pat(end+1) = struct('name', "8 cuff pads", 'ids', cuff8, 'w', eq(cuff8));
pat(end+1) = struct('name', "16 shank", 'ids', shank16, 'w', eq(shank16));
pat(end+1) = struct('name', "all 24 distal (wall)", 'ids', 25:48, 'w', eq(25:48));
pat(end+1) = struct('name', "nearest to F7, 4 stations", 'ids', nearest, 'w', eq(nearest));

%% Exposure per pattern ----------------------------------------------------------------------------
rows = cell(0, 4);
fig = figure('Visible', 'off', 'Position', [50 50 1500 820]);
tiledlayout(2, ceil(numel(pat) / 2), 'TileSpacing', 'compact');
for k = 1:numel(pat)
    tic;
    fx = fieldExposure(out, pat(k).ids, pat(k).w);
    nexttile;  hold on;
    for L = levels
        v = exposedVolume(fx, L, currents);
        vt = sum(v(:, 1:2), 2);                                % endoneurium + epineurium (the nerve)
        loglog(currents * 1e6, max(vt, 1e-9 * mm^3) / mm^3, 'DisplayName', sprintf('> %g V/cm', L / 100));
        for t = 1:numel(fx.tissue)                             % critical current per tissue
            n = find(v(:, t) >= Vcrit, 1, 'first');
            Ic = NaN;   if ~isempty(n), Ic = currents(n); end
            rows(end+1, :) = {pat(k).name, fx.tissue(t), L / 100, Ic * 1e6};   %#ok<SAGROW>
        end
        np = find(fx.periField_perA(1) * currents >= L, 1, 'first');   % perineurium: any node
        Ip = NaN;   if ~isempty(np), Ip = currents(np); end
        rows(end+1, :) = {pat(k).name, "perineurium (any node)", L / 100, Ip * 1e6};   %#ok<SAGROW>
    end
    set(gca, 'XScale', 'log', 'YScale', 'log');
    yline(Vcrit / mm^3, 'k:', 'HandleVisibility', 'off');
    xlabel('total current [\muA]');  ylabel('nerve volume above level [mm^3]');
    title(pat(k).name + sprintf(" (%d contacts%s)", numel(pat(k).ids), ...
          string(repmat(', mirrored', 1, any(fx.mirrored)))), 'FontSize', 9);
    ylim([1e-6 10]);  grid on;
    if k == 1, legend('Location', 'northwest', 'FontSize', 7); end
    fprintf('   %-28s done [%.0f s]\n', pat(k).name, toc);
end
sgtitle(sprintf('%s: volume of nerve tissue (endo + epineurium) where |E| exceeds each level; dotted = %g mm^3', ...
        b.patientId, Vcrit / mm^3), 'FontSize', 10);
exportgraphics(fig, fullfile(outDir, "field_exposure_" + b.patientId + ".png"), 'Resolution', 130);
close(fig);

T = cell2table(rows, 'VariableNames', {'pattern', 'tissue', 'level_Vpcm', 'criticalCurrent_uA'});
writetable(T, fullfile(outDir, "field_exposure_" + b.patientId + ".csv"));
fprintf('\nCritical current = total current at which %.0e mm^3 of the tissue exceeds the level (perineurium: any node):\n', ...
        Vcrit / mm^3);
disp(unstack(T, 'criticalCurrent_uA', 'level_Vpcm', 'VariableNamingRule', 'preserve'));
fprintf('Wrote %s\n', fullfile(outDir, "field_exposure_" + b.patientId + ".csv/png"));

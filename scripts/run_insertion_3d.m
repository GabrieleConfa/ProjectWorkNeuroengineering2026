function summary = run_insertion_3d(nTrials, shapesMode)
%RUN_INSERTION_3D Monte Carlo of the lead (probe) insertion on the patient of the main 3D model.
%   summary = run_insertion_3d()                     100 trials per station, circular fascicles
%   summary = run_insertion_3d(nTrials, shapesMode)  shapesMode = "circles" | "irregular"
%
% What it does
%   Adapted on 2026-11-26 from Marco's scripts/run_insertion_study.m (branch
%   cuff-oval-viewer-convergence; algorithm in docs/algorithms_and_plan.md §4), to run on OUR virtual
%   patient (makePatient + config/pw6_3d_params.m, the anatomy of main_pipeline_3d) through
%   src/insertion/insertionGeometry.m. For every one of the 8 stations, nTrials random insertions of
%   the 100 um PDMS lead are simulated in the station plane (2D), and two models of what happens when
%   the tip meets a fascicle are compared:
%     "force"  (Marco): grazing contact -> the lead slides around; head-on -> the fascicle is pushed
%              aside if there is room and the hand force beats its stiffness, otherwise the lead
%              pierces it if the force beats the puncture force, otherwise it slides. Random entry
%              angle and offset. The forces are NOT calibrated (config/insertionParams.m).
%     "probability" (tutor, meeting 2026-11-23, with the angle added by the user on 2026-11-26): at
%              the first contact with a fascicle, pierce with probability P = P_r(r) * cos(gamma)^m
%              (src/insertion/pierceProbability.m: small fascicles very unlikely, large ones likely,
%              less likely at a grazing angle; law and values [assumed], to ask the tutor), otherwise
%              bend around;
%              straight entry; after a bend the lead steers back to the entry line, so the exit
%              stays on the line of the entry.
%   It is a stand-alone study: it does NOT change the 3D geometry, the mesh or the solved contacts
%   (the 3D model still uses straight shanks, D-09/D-11; integrating the bent lead is the next step,
%   docs/design.md roadmap item 1).
%
% Outputs (results/insertion3d/)
%   fascicle_shapes.png                circles of the 3D model vs random irregular shapes
%   <model>_before_after.png           one station before and after an insertion
%   <model>_all_stations.png           one example insertion per station
%   <model>_statistics.png             outcome per station, fascicles crossed (+ force sensitivity)
%   pierce_probability.png             the assumed law P(r, angle), and P(pierce | contact) measured in the
%                                      Monte Carlo against fascicle radius and contact angle, both models
%   insertion_summary.csv, insertion_3d.mat
% Run: setup_paths; run_insertion_3d          Time: ~1-3 min (2D geometry only).

if nargin < 1 || isempty(nTrials), nTrials = []; end
if nargin < 2, shapesMode = "circles"; end
root = fileparts(fileparts(mfilename('fullpath')));   addpath(root);   setup_paths();
outDir = fullfile(root, 'results', 'insertion3d');   if ~isfolder(outDir), mkdir(outDir); end

%% 1. Patient of the main 3D model, in the format of the insertion code ---------------------------------
p = pw6_3d_params();
patient = makePatient(p.patient.id, p.patient.seed, p);
q = insertionGeometry(patient, p);
ip = insertionParams();   if ~isempty(nTrials), ip.nTrials = nTrials; end
th = (0:127)' * 2 * pi / 128;
circles = arrayfun(@(k) q.geometry.fascicles_m(k, 1:2) + q.geometry.fascicles_m(k, 3) * [cos(th) sin(th)], ...
                   (1:size(q.geometry.fascicles_m, 1))', 'UniformOutput', false);
irregular = randomFascicleShapes(q, ip, ip.seed);              % same centres and areas, random shapes
if shapesMode == "irregular", polys = irregular.polygons; else, polys = circles; end
shapes = struct('polygons', {polys}, 'nerve', irregular.nerve, ...
                'area_m2', cellfun(@(P) polyarea(P(:,1), P(:,2)), polys(:)), ...
                'roundness', ones(numel(polys), 1));
if shapesMode == "irregular", shapes.roundness = irregular.roundness; end
plotShapes(q, irregular, fullfile(outDir, 'fascicle_shapes.png'));
n = numel(polys);   rEq_mm = sqrt(shapes.area_m2 / pi) * 1e3;

%% 2. Stations from the provided contact layout (as in Marco's study) -------------------------------------
[st, stations] = insertionStations(q, shapes.nerve);
nS = numel(stations);

%% 3. Monte Carlo, both models ----------------------------------------------------------------------------
models = ["force" "probability"];
rows = table();   P = struct();
for m = models
    ipm = ip;   ipm.pierceModel = m;
    if m == "force"                                            % Marco's original settings of his force model
        ipm.entryAngleRange_deg = 10;   ipm.returnToLine = false;
    else                                                       % tutor's rule: straight entry, exit on its line
        ipm.entryAngleRange_deg = 0;   ipm.returnToLine = true;
    end
    counts = zeros(nS, 4);   pierced = zeros(1, n);   contacts = zeros(0, 3);   contacted = zeros(1, n);   ex = cell(nS, 1);   exT = zeros(nS, 1);
    padIn = 0;   padAll = 0;
    for iS = 1:nS
        R = simulateStation(ipm, q, polys, st{iS}, ipm.nTrials, ipm.seed + 17 * iS);
        counts(iS, :) = R.counts;   pierced = pierced + R.pierceCount;   contacted = contacted + R.contactCount;   contacts = [contacts; R.contacts]; %#ok<AGROW>
        ex{iS} = R.example;   exT(iS) = R.exampleTrial;   padIn = padIn + R.padIn;   padAll = padAll + R.padAll;
        fprintf('[%s] station %+d: no contact %d, slid only %d, displaced %d, pierced %d; completed %d/%d\n', ...
                m, stations(iS), R.counts, R.completed, ipm.nTrials);
    end
    frac = 100 * counts ./ sum(counts, 2);
    rows = [rows; table(repmat(m, nS, 1), stations(:), frac(:, 1), frac(:, 2), frac(:, 3), frac(:, 4), ...
                        'VariableNames', {'model', 'station', 'noContact_pct', 'slidOnly_pct', 'displaced_pct', 'pierced_pct'})]; %#ok<AGROW>
    P.(m) = struct('pierced', pierced, 'contacted', contacted, 'contacts', contacts, 'frac', frac, ...
                   'padsInFascicle_pct', 100 * padIn / padAll);
    plotBeforeAfter(shapes, ex{1}, stations(1), exT(1), m, fullfile(outDir, m + "_before_after.png"));
    plotAllStations(shapes, ex, stations, m, fullfile(outDir, m + "_all_stations.png"));
    sens = [];   forces = [];
    if m == "force"                                            % sensitivity to the uncalibrated hand force
        forces = [1 2 5 20 50] * 1e-3;   sens = zeros(numel(forces), 3);   nSens = min(ipm.nTrials, 25);
        for j = 1:numel(forces)
            ip2 = ipm;   ip2.pushForceMedian_N = forces(j);   c = zeros(1, 4);
            for iS = 1:nS, R = simulateStation(ip2, q, polys, st{iS}, nSens, ipm.seed + 101 * iS); c = c + R.counts; end
            sens(j, :) = 100 * c(2:4) / sum(c);
        end
    end
    plotStatistics(shapes, ipm, stations, frac / 100, pierced, forces, sens, m, fullfile(outDir, m + "_statistics.png"));
end
plotPierceProbability(rEq_mm, P, ip, fullfile(outDir, 'pierce_probability.png'));

summary = rows;
writetable(rows, fullfile(outDir, 'insertion_summary.csv'));
save(fullfile(outDir, 'insertion_3d.mat'), 'ip', 'shapes', 'shapesMode', 'stations', 'P', 'rows', 'rEq_mm');
disp(rows);
fprintf('Pads ending inside a fascicle: force model %.0f %%, probability model %.0f %%\n', ...
        P.force.padsInFascicle_pct, P.probability.padsInFascicle_pct);
fprintf('Wrote %s\n', outDir);
end


%% Local functions -------------------------------------------------------------------------------------------
function R = simulateStation(ip, q, polys, st, nTrials, seed)
% All trials of one station and their statistics (from Marco's run_insertion_study, plus the number
% of trials in which each fascicle was touched, for P(pierce | contact)).
n = numel(polys);   rs = RandStream('mt19937ar', 'Seed', seed);
rec = repmat(struct('deflect', false, 'displace', false, 'pierce', false, 'completed', false), nTrials, 1);
R = struct('pierceCount', zeros(1, n), 'contactCount', zeros(1, n), 'contacts', zeros(0, 3), ...
           'padIn', 0, 'padAll', 0, 'example', [], 'exampleTrial', 0);   % contacts: [r_m, angle_deg, pierced]
for t = 1:nTrials
    o = simulateInsertion(q, ip, polys, st, rs);
    rec(t).completed = o.completed;   rec(t).deflect = any(o.deflect);
    rec(t).displace = any(o.displace);   rec(t).pierce = any(o.pierce);
    R.pierceCount = R.pierceCount + double(o.pierce(:))';
    touched = o.pierce(:) | o.displace(:) | o.deflect(:);
    R.contactCount = R.contactCount + double(touched)';
    R.contacts = [R.contacts; o.rEq_m(touched), o.contactAngle_deg(touched), double(o.pierce(touched))];
    ok = ~isnan(o.padXY(:, 1));   R.padAll = R.padAll + nnz(ok);   R.padIn = R.padIn + nnz(o.padTissue(ok) > 0);
    if isempty(R.example) && o.completed && (any(o.displace) || any(o.pierce)), o.st = st; R.example = o; R.exampleTrial = t; end
end
if isempty(R.example), o.st = st; R.example = o; R.exampleTrial = nTrials; end
pi_ = [rec.pierce];   di = [rec.displace] & ~pi_;   de = [rec.deflect] & ~pi_ & ~di;   fr = ~[rec.deflect] & ~pi_ & ~di;
R.counts = [nnz(fr) nnz(de) nnz(di) nnz(pi_)];   R.completed = nnz([rec.completed]);
end

function plotShapes(q, shapes, file)
fig = figure(Visible="off", Color="w", Position=[40 40 1250 560]);   tl = tiledlayout(fig, 1, 2, TileSpacing="compact");
F = q.geometry.fascicles_m;   th = linspace(0, 2*pi, 200)';   nv = shapes.nerve;
ax = nexttile(tl);   hold(ax, 'on');   fill(ax, nv(:,1)*1e3, nv(:,2)*1e3, [0.97 0.91 0.77]);
for k = 1:size(F, 1)
    fill(ax, 1e3*(F(k,1) + F(k,3)*cos(th)), 1e3*(F(k,2) + F(k,4)*sin(th)), [0.96 0.75 0.68]);
    text(ax, F(k,1)*1e3, F(k,2)*1e3, "F" + k, HorizontalAlignment="center");
end
axis(ax, 'equal');   xlabel(ax, 'x (mm)');   ylabel(ax, 'y (mm)');   title(ax, 'Main 3D model: circular fascicles');
ax = nexttile(tl);   hold(ax, 'on');   fill(ax, nv(:,1)*1e3, nv(:,2)*1e3, [0.97 0.91 0.77]);
for k = 1:numel(shapes.polygons)
    Pk = shapes.polygons{k};   fill(ax, Pk(:,1)*1e3, Pk(:,2)*1e3, [0.96 0.75 0.68]);
    text(ax, mean(Pk(:,1))*1e3, mean(Pk(:,2))*1e3, sprintf('F%d', k), HorizontalAlignment="center");
end
axis(ax, 'equal');   xlabel(ax, 'x (mm)');   ylabel(ax, 'y (mm)');
title(ax, sprintf('Random irregular shapes (Marco): roundness %.2f \\pm %.2f (histology 0.68 \\pm 0.08)', ...
      mean(shapes.roundness), std(shapes.roundness)));
title(tl, 'Same centres and areas');   exportgraphics(fig, file, Resolution=150);   close(fig);
end

function drawSection(ax, shapes, P, origP)
nv = shapes.nerve;   hold(ax, 'on');   fill(ax, nv(:,1)*1e3, nv(:,2)*1e3, [0.97 0.91 0.77], EdgeColor=[0.3 0.3 0.3]);
for k = 1:numel(P)
    if ~isempty(origP), plot(ax, [origP{k}(:,1); origP{k}(1,1)]*1e3, [origP{k}(:,2); origP{k}(1,2)]*1e3, '--', Color=[0.5 0.5 0.5], LineWidth=0.8); end
    fill(ax, P{k}(:,1)*1e3, P{k}(:,2)*1e3, [0.96 0.75 0.68], EdgeColor=[0.65 0.40 0.35]);
    text(ax, mean(P{k}(:,1))*1e3, mean(P{k}(:,2))*1e3, "F" + k, HorizontalAlignment="center", FontSize=8);
end
lim = max(abs(nv), [], 1) * 1e3 * 1.05;
axis(ax, 'equal');   xlim(ax, [-lim(1) lim(1)]);   ylim(ax, [-lim(2) lim(2)]);   xlabel(ax, 'x (mm)');   ylabel(ax, 'y (mm)');
end

function drawLead(ax, o, st)
plot(ax, o.path(:,1)*1e3, o.path(:,2)*1e3, '-', Color=[0.12 0.36 0.94], LineWidth=1.8);
for s = 1:size(o.pierceSeg, 1)
    idx = o.pierceSeg(s,1):o.pierceSeg(s,2);   plot(ax, o.path(idx,1)*1e3, o.path(idx,2)*1e3, '-', Color=[0.85 0.1 0.1], LineWidth=3);
end
for k = 1:numel(st.padS)
    if isnan(o.padXY(k,1)), continue; end
    e = o.padXY(k,:) + [-1; 1] * st.padHalf * o.padTan(k,:);   plot(ax, e(:,1)*1e3, e(:,2)*1e3, '-', Color=[0.62 0.13 0.85], LineWidth=4);
end
end

function plotBeforeAfter(shapes, o, sec, trial, model, file)
st = o.st;   fig = figure(Visible="off", Color="w", Position=[40 40 1300 620]);   tl = tiledlayout(fig, 1, 2, TileSpacing="compact");
ax = nexttile(tl);   drawSection(ax, shapes, shapes.polygons, []);   hold(ax, 'on');
a = [st.entryS; st.tipS] * st.u;   plot(ax, a(:,1)*1e3, a(:,2)*1e3, '-', Color=[0.12 0.36 0.94], LineWidth=1.8);
for k = 1:numel(st.padS), c = st.padS(k) * st.u; e = c + [-1; 1] * st.padHalf * st.u; plot(ax, e(:,1)*1e3, e(:,2)*1e3, '-', Color=[0.62 0.13 0.85], LineWidth=4); end
title(ax, sprintf('Before insertion (station %+d): nominal straight lead', sec));
ax = nexttile(tl);   drawSection(ax, shapes, o.polys, shapes.polygons);   drawLead(ax, o, st);
ev = "";
if any(o.pierce), ev = ev + " pierced: F" + strjoin(string(find(o.pierce)'), ",F") + ";"; end
if any(o.displace), ev = ev + " displaced: F" + strjoin(string(find(o.displace)'), ",F") + ";"; end
if any(o.deflect), ev = ev + " slid around: F" + strjoin(string(find(o.deflect)'), ",F") + ";"; end
title(ax, sprintf('After insertion (trial %d):%s', trial, ev));
title(tl, sprintf('Model "%s": fascicles before and after one random insertion (dashed: original position; red: pierced stretch; purple: pads)', model));
exportgraphics(fig, file, Resolution=150);   close(fig);
end

function plotAllStations(shapes, ex, stations, model, file)
fig = figure(Visible="off", Color="w", Position=[40 40 1500 640]);   tl = tiledlayout(fig, 2, 4, TileSpacing="compact");
for iS = 1:numel(stations)
    ax = nexttile(tl);   o = ex{iS};   drawSection(ax, shapes, o.polys, shapes.polygons);   drawLead(ax, o, o.st);
    ev = "";   if any(o.pierce), ev = ev + "P"; end, if any(o.displace), ev = ev + "D"; end, if any(o.deflect), ev = ev + "S"; end
    title(ax, sprintf('Station %+d  [%s]', stations(iS), ev), FontSize=10);
end
title(tl, sprintf('Model "%s": one random insertion per station (P pierce, D displace, S slide; dashed = original position)', model));
exportgraphics(fig, file, Resolution=140);   close(fig);
end

function plotStatistics(shapes, ip, stations, frac, pierceCount, forces, sens, model, file)
nPan = 2 + ~isempty(forces);
fig = figure(Visible="off", Color="w", Position=[40 40 570*nPan 520]);   tl = tiledlayout(fig, 1, nPan, TileSpacing="compact");
ax = nexttile(tl);   b = bar(ax, frac*100, 'stacked');   cols = [0.75 0.85 0.95; 0.35 0.65 0.85; 0.95 0.75 0.25; 0.85 0.15 0.15];
for k = 1:4, b(k).FaceColor = cols(k,:); end
set(ax, XTickLabel=compose('%+d', stations));   xlabel(ax, 'Station');   ylabel(ax, 'Insertions (%)');
legend(ax, {'no contact', 'slid around only', 'fascicle displaced', 'fascicle pierced'}, Location='southoutside', Orientation='horizontal');
grid(ax, 'on');   title(ax, sprintf('Model "%s": outcome of %d insertions per station', model, ip.nTrials));
ax = nexttile(tl);   a = shapes.area_m2*1e6;   tot = pierceCount / (ip.nTrials*numel(stations)) * 100;
[~, o] = sort(a, 'descend');   bar(ax, tot(o));
set(ax, XTick=1:numel(a), XTickLabel=arrayfun(@(k) sprintf('F%d\n%.2f', o(k), a(o(k))), 1:numel(a), 'UniformOutput', false));
grid(ax, 'on');   xlabel(ax, 'Fascicle and area (mm^2)');   ylabel(ax, 'Trials in which it was pierced (%)');   title(ax, 'Which fascicles are crossed');
if ~isempty(forces)
    ax = nexttile(tl);   hold(ax, 'on');   grid(ax, 'on');
    plot(ax, forces*1e3, sens(:,3), '-o', LineWidth=1.8, DisplayName='pierced');
    plot(ax, forces*1e3, sens(:,2), '-s', LineWidth=1.8, DisplayName='displaced');
    plot(ax, forces*1e3, sens(:,1), '-d', LineWidth=1.8, DisplayName='slid only');
    set(ax, XScale='log');   xlabel(ax, 'Median hand force (mN)');   ylabel(ax, 'Trials (%)');   legend(ax, Location='best');
    title(ax, 'Sensitivity to the uncalibrated hand force');
end
exportgraphics(fig, file, Resolution=150);   close(fig);
end

function plotPierceProbability(rEq_mm, P, ip, file)
% (a) the assumed law P(r, angle) of the probability model, with the radii of the patient's fascicles;
% (b) P(pierce | contact) per fascicle against its radius, measured in the Monte Carlo, both models;
% (c) P(pierce | contact) against the angle of the first contact, measured, both models.
fig = figure(Visible="off", Color="w", Position=[40 40 1700 520]);   tl = tiledlayout(fig, 1, 3, TileSpacing="compact");
% (a) the law
ax = nexttile(tl);
r = linspace(0.1, 0.8, 141);   g = linspace(0, 90, 91);   [RR, GG] = meshgrid(r, g);
PP = 100 * pierceProbability(RR * 1e-3, GG, ip);
imagesc(ax, r, g, PP);   set(ax, 'YDir', 'normal');   hold(ax, 'on');
colormap(ax, parula);   cb = colorbar(ax);   cb.Label.String = 'P(pierce) (%)';
contour(ax, r, g, PP, [10 50 90], 'k', 'ShowText', 'on');
for k = 1:numel(rEq_mm)
    xline(ax, rEq_mm(k), 'w:');
    text(ax, rEq_mm(k), 87, "F" + k, 'Color', 'w', 'FontSize', 7, 'HorizontalAlignment', 'center');
end
xlabel(ax, 'fascicle radius sqrt(A/\pi) (mm)');   ylabel(ax, 'contact angle to the normal (deg; 0 = head-on)');
title(ax, sprintf('Assumed law: P = 1/(1+(r_{50}/r)^n) cos^m(\\gamma),  r_{50} = %.2f mm, n = %g, m = %g', ...
      ip.pierceR50_m * 1e3, ip.pierceSteepness, ip.pierceAngleExponent), 'FontSize', 9);
% (b) measured, against the radius
ax = nexttile(tl);   hold(ax, 'on');   grid(ax, 'on');
mk = struct('force', 'o', 'probability', 's');   col = struct('force', [0.12 0.36 0.94], 'probability', [0.85 0.15 0.15]);
for m = ["force" "probability"]
    c = P.(m).contacted;   ok = c > 0;
    scatter(ax, rEq_mm(ok), 100 * P.(m).pierced(ok) ./ c(ok), 20 + 2 * c(ok), col.(m), mk.(m), 'filled', ...
            'MarkerFaceAlpha', 0.6, 'DisplayName', sprintf('measured, model "%s" (size ~ contacts)', m));
end
rr = linspace(0.1, max(rEq_mm) * 1.1, 200);
for gdeg = [0 45]
    plot(ax, rr, 100 * pierceProbability(rr * 1e-3, gdeg, ip), 'k--', 'LineWidth', 1 + (gdeg == 0), ...
         'DisplayName', sprintf('law at %d deg', gdeg));
end
xlabel(ax, 'fascicle radius (mm)');   ylabel(ax, 'P(pierce | contact) (%)');   ylim(ax, [-5 105]);
legend(ax, 'Location', 'northwest', 'FontSize', 7);
title(ax, 'Measured per fascicle (all contact angles together)', 'FontSize', 9);
% (c) measured, against the contact angle
ax = nexttile(tl);   hold(ax, 'on');   grid(ax, 'on');
edges = 0:15:90;   mid = edges(1:end-1) + 7.5;
for m = ["force" "probability"]
    C = P.(m).contacts;   if isempty(C), continue; end
    b = discretize(C(:, 2), edges);   pr = nan(size(mid));
    for j = 1:numel(mid), sel = b == j; if nnz(sel) >= 5, pr(j) = 100 * mean(C(sel, 3)); end, end
    plot(ax, mid, pr, ['-' mk.(m)], 'Color', col.(m), 'LineWidth', 1.6, 'MarkerFaceColor', col.(m), ...
         'DisplayName', sprintf('measured, model "%s"', m));
end
xlabel(ax, 'contact angle to the normal (deg; 0 = head-on)');   ylabel(ax, 'P(pierce | contact) (%)');   ylim(ax, [-5 105]);
legend(ax, 'Location', 'northeast', 'FontSize', 7);
title(ax, 'Measured against the angle of the first contact (all fascicles; bins with >= 5 contacts)', 'FontSize', 9);
title(tl, {'Probability that the lead pierces a fascicle it touches', ...
           'indicative only: the forces of the force model and the law of the probability model are not calibrated'});
exportgraphics(fig, file, Resolution=150);   close(fig);
end

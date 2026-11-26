function summary = run_insertion_study(nTrials) % Monte Carlo of random lead insertions with irregular random fascicles.
%% 0. What this study does
% 1) draws random irregular fascicle outlines (roundness statistics from histology),
% 2) simulates nTrials random insertions per station (deflect / displace / pierce, see simulateInsertion),
% 3) draws the section BEFORE and AFTER one insertion, summary statistics and the pad tissues.
% It is a stand-alone scenario study: it does not change the production geometry, mesh or viewer.
% All parameters (and their status: GIVEN / LIT / ASSUMED) are in config/insertionParams.m.
root = fileparts(fileparts(mfilename('fullpath'))); addpath(root); setup_paths(); % Resolve all sources relative to this file.
outDir = fullfile(root,'results','insertion'); if ~isfolder(outDir), mkdir(outDir); end % Output folder.
q = quicklook3DParams(); ip = insertionParams(); if nargin>=1, ip.nTrials = nTrials; end % Anatomy and parameters.
%% 1. Random irregular fascicles
shapes = randomFascicleShapes(q,ip,ip.seed); polys = shapes.polygons; n = numel(polys); % Same area and centre as the circular layout.
plotShapes(q,shapes,fullfile(outDir,'fascicle_shapes.png')); % Circles versus random shapes.
%% 2. Stations from the supplied device file
raw = loadActiveSites(q.electrode.sitesFile,q.geometry.radius_m); stations = unique(raw.section)'; % Eight axial stations.
nS = numel(stations); st = cell(nS,1); % Station descriptions.
for iS = 1:nS
    sec = stations(iS); rows = find(raw.section==sec & raw.type=="intraneural"); sites = raw(rows,:); % The four internal pads.
    [~,i] = max(sites.r_m); u = sites.centre_m(i,1:2)/sites.r_m(i); % Nominal insertion axis (as in equilibriumStrip).
    if ismember(sec,q.electrode.reverseEntrySections), u = -u; end % Same entry side as the approved lead trajectories.
    padS = sites.centre_m(:,1:2)*u'; tipS = nerveBoundaryRadius(q,u); % Nominal pad positions; the lead ends at the opposite wall, where it is anchored.
    entryS = -nerveBoundaryRadius(q,-u); % Nominal entry wall position (negative side).
    st{iS} = struct('u',u,'padS',padS,'tipS',tipS,'entryS',entryS,'nerve',shapes.nerve,'section',sec,'padHalf',sites.sideXY_m(1)/2); % Station description.
end
counts = zeros(nS,4); pierceCount = zeros(nS,n); padIn = zeros(nS,1); padAll = zeros(nS,1); example = cell(nS,1); exampleTrial = zeros(nS,1); % Statistics.
for iS = 1:nS % Main Monte Carlo with the default parameters.
    R = simulateStation(q,ip,polys,st{iS},ip.nTrials,ip.seed+17*iS); % All trials of this station.
    counts(iS,:) = R.counts; pierceCount(iS,:) = R.pierceCount; padIn(iS) = R.padIn; padAll(iS) = R.padAll; example{iS} = R.example; exampleTrial(iS) = R.exampleTrial; % Store.
    fprintf('Station %+d: free %d, deflected %d, displaced %d, pierced %d, completed %d/%d\n',stations(iS),R.counts,R.completed,ip.nTrials); % Progress.
end
%% 2b. Sensitivity to the (uncalibrated) hand force
forces = [1 2 5 20 50]*1e-3; sens = zeros(numel(forces),3); nSens = min(ip.nTrials,25); % Median hand forces tested (N), pooled over stations.
for j = 1:numel(forces)
    ip2 = ip; ip2.pushForceMedian_N = forces(j); c = zeros(1,4); % Same parameters except the hand force.
    for iS = 1:nS, R = simulateStation(q,ip2,polys,st{iS},nSens,ip.seed+101*iS); c = c+R.counts; end % Pool the stations.
    sens(j,:) = 100*[c(2) c(3) c(4)]/sum(c); % Percent of trials: slid only, displaced, pierced.
    fprintf('Hand force median %.0f mN: slid %.0f %%, displaced %.0f %%, pierced %.0f %%\n',forces(j)*1e3,sens(j,:)); % Progress.
end
writetable(array2table([forces(:)*1e3 sens],'VariableNames',{'handForceMedian_mN','slid_pct','displaced_pct','pierced_pct'}),fullfile(outDir,'insertion_sensitivity.csv')); % Save.
%% 3. Tables
frac = counts./sum(counts,2); % Fractions of trials.
summary = array2table([stations(:) frac*100 padIn./padAll*100],'VariableNames',{'station','free_pct','deflectedOnly_pct','displaced_pct','pierced_pct','padsInsideFascicle_pct'}); % One row per station.
writetable(summary,fullfile(outDir,'insertion_summary.csv')); disp(summary); % Save and print.
save(fullfile(outDir,'insertion_study.mat'),'ip','shapes','counts','pierceCount','example','exampleTrial','stations','forces','sens'); % Keep everything.
%% 4. Figures
plotBeforeAfter(q,shapes,example{1},stations(1),exampleTrial(1),fullfile(outDir,'insertion_before_after.png')); % Detailed station -4 (first station).
plotAllStations(q,shapes,example,stations,fullfile(outDir,'insertion_all_stations.png')); % Overview of all stations.
plotStatistics(shapes,ip,stations,frac,pierceCount,forces,sens,fullfile(outDir,'insertion_statistics.png')); % Outcomes and fascicle size.
end
%% Local functions
function R = simulateStation(q,ip,polys,st,nTrials,seed) % Run all trials of one station and collect the statistics.
n = numel(polys); rs = RandStream('mt19937ar','Seed',seed); % Reproducible stream.
rec = repmat(struct('deflect',false,'displace',false,'pierce',false,'completed',false),nTrials,1); % Trial flags.
R = struct('pierceCount',zeros(1,n),'padIn',0,'padAll',0,'example',[],'exampleTrial',0); % Outputs.
for t = 1:nTrials
    o = simulateInsertion(q,ip,polys,st,rs); % One random insertion.
    rec(t).completed = o.completed; rec(t).deflect = any(o.deflect); rec(t).displace = any(o.displace); rec(t).pierce = any(o.pierce); % Flags.
    R.pierceCount = R.pierceCount+double(o.pierce(:))'; % Which fascicles were crossed.
    ok = ~isnan(o.padXY(:,1)); R.padAll = R.padAll+nnz(ok); R.padIn = R.padIn+nnz(o.padTissue(ok)>0); % Pads that ended inside a fascicle.
    if isempty(R.example) && o.completed && (any(o.displace) || any(o.pierce)), o.st = st; R.example = o; R.exampleTrial = t; end % First illustrative trial.
end
if isempty(R.example), o.st = st; R.example = o; R.exampleTrial = nTrials; end % Fall back to the last trial.
pi_ = [rec.pierce]; di = [rec.displace]&~pi_; de = [rec.deflect]&~pi_&~di; fr = ~[rec.deflect]&~pi_&~di; % Exclusive classes per trial.
R.counts = [nnz(fr) nnz(de) nnz(di) nnz(pi_)]; R.completed = nnz([rec.completed]); % Free, slid only, displaced (no pierce), pierced.
end
%% Local functions: figures
function plotShapes(q,shapes,file) % Original circles and random irregular shapes side by side.
fig = figure(Visible="off",Color="w",Position=[40 40 1250 560],Theme="light"); tl = tiledlayout(fig,1,2,TileSpacing="compact"); % Two panels.
F = q.geometry.fascicles_m; th = linspace(0,2*pi,200)'; nv = shapes.nerve; % Circles and nerve outline.
ax = nexttile(tl); hold(ax,'on'); fill(ax,nv(:,1)*1e3,nv(:,2)*1e3,[0.97 0.91 0.77]); % Nerve.
for k = 1:size(F,1), fill(ax,1e3*(F(k,1)+F(k,3)*cos(th)),1e3*(F(k,2)+F(k,4)*sin(th)),[0.96 0.75 0.68]); text(ax,F(k,1)*1e3,F(k,2)*1e3,"F"+k,HorizontalAlignment="center"); end % Circular fascicles.
axis(ax,'equal'); xlabel(ax,'x (mm)'); ylabel(ax,'y (mm)'); title(ax,'Previous layout: circular fascicles'); % Labels.
ax = nexttile(tl); hold(ax,'on'); fill(ax,nv(:,1)*1e3,nv(:,2)*1e3,[0.97 0.91 0.77]); % Nerve.
for k = 1:numel(shapes.polygons), P = shapes.polygons{k}; fill(ax,P(:,1)*1e3,P(:,2)*1e3,[0.96 0.75 0.68]); text(ax,mean(P(:,1))*1e3,mean(P(:,2))*1e3,sprintf('F%d',k),HorizontalAlignment="center"); end % Random shapes.
axis(ax,'equal'); xlabel(ax,'x (mm)'); ylabel(ax,'y (mm)'); % Labels.
title(ax,sprintf('Random irregular fascicles: roundness %.2f \\pm %.2f (histology 0.68 \\pm 0.08)',mean(shapes.roundness),std(shapes.roundness))); % Statistics.
title(tl,'Same centres and areas, random shapes'); exportgraphics(fig,file,Resolution=170); close(fig); % Save.
end
function drawSection(ax,shapes,P,origP) % Nerve, original outlines (dashed) and current fascicles.
nv = shapes.nerve; hold(ax,'on'); fill(ax,nv(:,1)*1e3,nv(:,2)*1e3,[0.97 0.91 0.77],EdgeColor=[0.3 0.3 0.3]); % Nerve.
for k = 1:numel(P) % Fascicles.
    if ~isempty(origP), plot(ax,[origP{k}(:,1);origP{k}(1,1)]*1e3,[origP{k}(:,2);origP{k}(1,2)]*1e3,'--',Color=[0.5 0.5 0.5],LineWidth=0.8); end % Original position.
    fill(ax,P{k}(:,1)*1e3,P{k}(:,2)*1e3,[0.96 0.75 0.68],EdgeColor=[0.65 0.40 0.35]); text(ax,mean(P{k}(:,1))*1e3,mean(P{k}(:,2))*1e3,"F"+k,HorizontalAlignment="center",FontSize=8); % Moved outline.
end
axis(ax,'equal'); xlim(ax,[-2.4 2.4]); ylim(ax,[-1.6 1.6]); xlabel(ax,'x (mm)'); ylabel(ax,'y (mm)'); % Common limits.
end
function drawLead(ax,o,st) % Lead path, pierced segments and pads.
plot(ax,o.path(:,1)*1e3,o.path(:,2)*1e3,'-',Color=[0.12 0.36 0.94],LineWidth=1.8); % Lead centreline.
for s = 1:size(o.pierceSeg,1), idx = o.pierceSeg(s,1):o.pierceSeg(s,2); plot(ax,o.path(idx,1)*1e3,o.path(idx,2)*1e3,'-',Color=[0.85 0.1 0.1],LineWidth=3); end % Pierced stretches.
for k = 1:numel(st.padS) % Pads as short segments along the lead.
    if isnan(o.padXY(k,1)), continue; end
    e = o.padXY(k,:)+[-1;1]*st.padHalf*o.padTan(k,:); plot(ax,e(:,1)*1e3,e(:,2)*1e3,'-',Color=[0.62 0.13 0.85],LineWidth=4); % Active face.
end
end
function plotBeforeAfter(q,shapes,o,sec,trial,file) % The same section before and after one insertion.
st = o.st; fig = figure(Visible="off",Color="w",Position=[40 40 1300 620],Theme="light"); tl = tiledlayout(fig,1,2,TileSpacing="compact"); % Two panels.
ax = nexttile(tl); drawSection(ax,shapes,shapes.polygons,[]); hold(ax,'on'); % Before.
a = [st.entryS;st.tipS]*st.u; plot(ax,a(:,1)*1e3,a(:,2)*1e3,'-',Color=[0.12 0.36 0.94],LineWidth=1.8); % Nominal straight lead.
for k = 1:numel(st.padS), c = st.padS(k)*st.u; tg = [-st.u(2) st.u(1)]*0; e = c+[-1;1]*st.padHalf*st.u; plot(ax,e(:,1)*1e3,e(:,2)*1e3,'-',Color=[0.62 0.13 0.85],LineWidth=4); end % Nominal pads.
title(ax,sprintf('Before insertion (station %+d): nominal straight lead',sec)); % Title.
ax = nexttile(tl); drawSection(ax,shapes,o.polys,shapes.polygons); drawLead(ax,o,st); % After.
ev = ""; if any(o.pierce), ev = ev+" pierced: F"+strjoin(string(find(o.pierce)'),",F")+";"; end % Events in the title.
if any(o.displace), ev = ev+" displaced: F"+strjoin(string(find(o.displace)'),",F")+";"; end
if any(o.deflect), ev = ev+" slid around: F"+strjoin(string(find(o.deflect)'),",F")+";"; end
title(ax,sprintf('After insertion (trial %d):%s',trial,ev)); % Title.
title(tl,'Fascicles before and after one random lead insertion (dashed: original position; red: pierced stretch; purple: pads)'); exportgraphics(fig,file,Resolution=170); close(fig); % Save.
end
function plotAllStations(q,shapes,ex,stations,file) % One example insertion per station.
fig = figure(Visible="off",Color="w",Position=[40 40 1500 640],Theme="light"); tl = tiledlayout(fig,2,4,TileSpacing="compact"); % 2 x 4 grid.
for iS = 1:numel(stations)
    ax = nexttile(tl); o = ex{iS}; drawSection(ax,shapes,o.polys,shapes.polygons); drawLead(ax,o,o.st); % Section after insertion.
    ev = ""; if any(o.pierce), ev = ev+"P"; end, if any(o.displace), ev = ev+"D"; end, if any(o.deflect), ev = ev+"S"; end % P pierce, D displace, S slide.
    title(ax,sprintf('Station %+d  [%s]',stations(iS),ev),FontSize=10); % Title.
end
title(tl,'One random insertion per station (P pierce, D displace, S slide; dashed = original fascicle position)'); exportgraphics(fig,file,Resolution=160); close(fig); % Save.
end
function plotStatistics(shapes,ip,stations,frac,pierceCount,forces,sens,file) % Outcome fractions and which fascicles are crossed.
fig = figure(Visible="off",Color="w",Position=[40 40 1700 520],Theme="light"); tl = tiledlayout(fig,1,3,TileSpacing="compact"); % Three panels.
ax = nexttile(tl); b = bar(ax,frac*100,'stacked'); cols = [0.75 0.85 0.95; 0.35 0.65 0.85; 0.95 0.75 0.25; 0.85 0.15 0.15]; % Four classes.
for k = 1:4, b(k).FaceColor = cols(k,:); end
set(ax,XTickLabel=compose('%+d',stations)); xlabel(ax,'Station'); ylabel(ax,'Insertions (%)'); legend(ax,{'no contact','slid around only','fascicle displaced','fascicle pierced'},Location='southoutside',Orientation='horizontal'); grid(ax,'on'); % Labels.
title(ax,sprintf('Outcome of %d random insertions per station',ip.nTrials)); % Title.
ax = nexttile(tl); a = shapes.area_m2*1e6; tot = sum(pierceCount,1)/(ip.nTrials*numel(stations))*100; % Pierce rate per fascicle.
[~,o] = sort(a,'descend'); bar(ax,tot(o)); set(ax,XTick=1:numel(a),XTickLabel=arrayfun(@(k)sprintf('F%d\n%.2f',o(k),a(o(k))),1:numel(a),'UniformOutput',false)); grid(ax,'on'); % Sorted by area.
xlabel(ax,'Fascicle and area (mm^2)'); ylabel(ax,'Trials in which it was pierced (%)'); title(ax,'Which fascicles are crossed'); % Labels.
ax = nexttile(tl); hold(ax,'on'); grid(ax,'on'); plot(ax,forces*1e3,sens(:,3),'-o',LineWidth=1.8,DisplayName='pierced'); plot(ax,forces*1e3,sens(:,2),'-s',LineWidth=1.8,DisplayName='displaced'); plot(ax,forces*1e3,sens(:,1),'-d',LineWidth=1.8,DisplayName='slid only'); % Sensitivity.
set(ax,XScale='log'); xlabel(ax,'Median hand force (mN)'); ylabel(ax,'Trials (%)'); legend(ax,Location='best'); title(ax,'Sensitivity to the uncalibrated hand force'); % Labels.
exportgraphics(fig,file,Resolution=170); close(fig); % Save.
end

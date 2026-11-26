function figures = plotStationAtlas(model) % Draw every electrode section and one close view of strip deflection.
%% 1. Prepare the current anatomical cross-section
p = model.parameters; % Read the editable nerve and fascicle dimensions.
mm = 1e3; % Convert model coordinates from metres to millimetres.
fasc = p.geometry.fascicles_m; % Read fascicle centres, axes and angles.
theta = linspace(0,2*pi,241); % Sample smooth anatomical boundaries.
figures.atlas = figure(Name="Eight electrode sections",Color="w",Theme="light",Position=[60 70 1420 800]); % Create a page for all eight levels.
layout = tiledlayout(figures.atlas,2,4,TileSpacing="compact",Padding="compact"); % Reserve one panel per station.
%% 2. Plot all eight measured axial stations
for k = 1:numel(model.strips) % Walk through every axial section in the tutor's file.
    strip = model.strips(k); % Read the planned strip and finite active faces.
    ax = nexttile(layout); % Select this station's panel.
    hold(ax,'on'); % Layer tissue regions and contacts in the same section.
    z_m = strip.path_m(1,3); % Read the physical position of the station along the nerve.
    [~,iz] = min(abs(model.nerve.z_m(:,1)-z_m)); % Find the sampled nerve boundary at this station.
    x = mm*model.nerve.x_m(iz,:); % Read the locally flattened outer contour.
    y = mm*model.nerve.y_m(iz,:); % Read its transverse y positions.
    epi = fill(ax,x,y,[0.97 0.91 0.77],EdgeColor=[0.35 0.29 0.22],LineWidth=1.5); % Show epineurium and its outer edge.
    plot(ax,x,y,Color=[0.13 0.38 0.92],LineWidth=2); % Draw the cuff wrapped around this station.
    fp = stationFasc(strip,p); h = drawFascicles(ax,fp,p.geometry.tPeri_m,mm,true); peri = h.peri; endo = h.endo; % Fascicles of this station (random shapes after the insertion).
    [segments_m,straight_m] = visibleInsertionSegments(strip); % Read the single-entry carrier.
    path = mm*segments_m(:,1:2); % Convert the connected lead to millimetres.
    reference = plot(ax,mm*straight_m(:,1),mm*straight_m(:,2),'--',Color=[0.43 0.43 0.43],LineWidth=1); % Show the undeviated insertion axis.
    device = plot(ax,path(:,1),path(:,2),Color=[0.13 0.38 0.92],LineWidth=2); % Show only separate insertion segments.
    for j = 1:size(strip.padFaces_m,3) % Show all four internal active faces at this level.
        face = mm*strip.padFaces_m(:,:,j); % Read the finite pad vertices.
        pad = plot(ax,face(1:2,1),face(1:2,2),Color=[0.57 0.14 0.83],LineWidth=5); % Project the active pad width.
    end % Finish the internal pads.
    entry = plot(ax,path(1,1),path(1,2),'ko',MarkerFaceColor='w',MarkerSize=5); % Mark the entry so handedness can be read in the insertion frame.
    cuff = model.nerve.cuffFaces_m; % Read all external pad faces.
    for j = 1:size(cuff,3) % Find the two cuff contacts at this station.
        face = mm*cuff(:,:,j); % Convert cuff vertices to millimetres.
        if abs(mean(face(:,3))-mm*z_m)<1e-5 % Select the current axial level.
            plot(ax,face(1:2,1),face(1:2,2),Color=[0.57 0.14 0.83],LineWidth=6); % Show the opposing cuff pad.
        end % Skip contacts at other levels.
    end % Finish locating cuff contacts.
    axis(ax,'equal'); % Preserve anatomical proportions.
    xlim(ax,mm*max(p.geometry.nerveSemiAxes_m)*[-1.12 1.12]); % Use the same x scale in all panels.
    ylim(ax,mm*max(p.geometry.nerveSemiAxes_m)*[-1.12 1.12]); % Use the same y scale in all panels.
    title(ax,sprintf('Section %+d  |  z = %.1f mm',strip.section,mm*z_m)); % Give the physical location and station ID.
    xlabel(ax,'x [mm]'); % Label transverse x position.
    ylabel(ax,'y [mm]'); % Label transverse y position.
end % Finish the atlas of all stations.
title(layout,'Median nerve cross-sections'); % Keep the shared title independent of the colour key.
key = legend(ax,[epi peri endo device pad reference entry], ... % Use actual drawn objects so each symbol matches the atlas.
    {'Epineurium','Perineurium','Endoneurium','Device (cuff and lead)','Electrode pads','Insertion axis','Entry point'}, ... % Identify the wrapping device and both sets of electrodes.
    Orientation="horizontal",NumColumns=4,Box="off"); % Arrange the seven entries in two rows.
key.Layout.Tile='south'; % Reserve space below all eight sections for the common legend.
%% 3. Highlight one real detour around the central fascicle
figures.deviation = figure(Name="Strip deflection detail",Color="w",Theme="light",Position=[130 120 870 760]); % Create a large explanation view.
ax = axes(figures.deviation); % Create axes for a single section.
hold(ax,'on'); % Overlay anatomy, nominal path and routed path.
strip = model.strips(1); % Select section -4 for the close-up.
[~,iz] = min(abs(model.nerve.z_m(:,1)-strip.path_m(1,3))); % Select its flattened cross-section.
fill(ax,mm*model.nerve.x_m(iz,:),mm*model.nerve.y_m(iz,:),[0.97 0.91 0.77], ...
    EdgeColor=[0.35 0.29 0.22],LineWidth=1.6); % Show the epineurium and outer edge.
drawFascicles(ax,stationFasc(strip,p),p.geometry.tPeri_m,mm,true); % Fascicles of this station.
[segments_m,straight_m] = visibleInsertionSegments(strip); % Get the routed and straight single-entry paths.
path = mm*segments_m(:,1:2); % Read only the visible candidate segments in millimetres.
straight = mm*straight_m(:,1:2); % Read their straight reference without a through-nerve line.
routedLine = plot(ax,path(:,1),path(:,2),Color=[0.05 0.40 0.62],LineWidth=3,DisplayName="Candidate detour"); % Show the exploratory separated routes.
straightLine = plot(ax,straight(:,1),straight(:,2),'--', ...
    Color=[0.78 0.19 0.17],LineWidth=1.8,DisplayName="Straight reference"); % Compare only the partial insertion lengths.
for j = 1:size(model.nerve.cuffFaces_m,3) % Show the external contacts at this selected section.
    cuffFace = mm*model.nerve.cuffFaces_m(:,:,j); % Convert this cuff face to millimetres.
    if abs(mean(cuffFace(:,3))-mm*strip.path_m(1,3))<1e-5 % Select this section's opposing cuff pair.
        plot(ax,cuffFace(1:2,1),cuffFace(1:2,2),Color=[0.57 0.14 0.83],LineWidth=6); % Mark the rigid contact surface.
    end % Skip the other axial stations.
end % Finish displaying the external contacts.
for j = 1:size(strip.padFaces_m,3) % Show each internal contact on its routed position.
    face = mm*strip.padFaces_m(:,:,j); % Read this finite active face.
    plot(ax,face(1:2,1),face(1:2,2),Color=[0.57 0.14 0.83],LineWidth=6); % Highlight the contact width.
end % Finish drawing internal active sites.
axis(ax,'equal'); % Keep the comparison anatomically proportionate.
lim = mm*max(p.geometry.nerveSemiAxes_m)*1.12; % Frame the full nerve section.
xlim(ax,[-lim lim]); ylim(ax,[-lim lim]); % Use equal transverse viewing limits.
xlabel(ax,'x [mm]'); ylabel(ax,'y [mm]'); % State the position units.
title(ax,sprintf('Section %+d | candidate detour %.3f mm | single-entry lead',strip.section,mm*abs(strip.offset_m))); % Mark the provisional geometric detour.
legend(ax,[routedLine,straightLine],{'Routed lead','Straight reference'},'Location','southoutside'); % Explain the two schematic comparisons.
end % Return both anatomy figures.
function xy = ellipsePoints(f,peri_m,theta) % Sample one rotated fascicle and optional sheath.
a = f(3)+peri_m; b = f(4)+peri_m; % Read both semi-axes with the outer perineurial thickness.
c = cosd(f(5)); s = sind(f(5)); % Rotate the fascicle in the transverse plane.
xy = [f(1)+a*cos(theta(:))*c-b*sin(theta(:))*s, ...
      f(2)+a*cos(theta(:))*s+b*sin(theta(:))*c]; % Return boundary coordinates in metres.
end % Finish the ellipse sampler.

function fp = stationFasc(strip,p) % Outlines of the station: from the strip (insertion mode) or from the editable ellipses.
if isfield(strip,'fascicles') && ~isempty(strip.fascicles), fp = strip.fascicles; else, fp = stationFascicles(p,strip.section); end
end

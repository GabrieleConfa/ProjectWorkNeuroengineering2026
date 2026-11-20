function figures = plotStationAtlas(model) % Draw every electrode section and one close view of strip deflection.
%% 1. Prepare the current anatomical cross-section
p = model.parameters; % Read the editable nerve and fascicle dimensions.
mm = 1e3; % Convert model coordinates from metres to millimetres.
fasc = p.geometry.fascicles_m; % Read fascicle centres, axes and angles.
theta = linspace(0,2*pi,241); % Sample smooth anatomical boundaries.
colors = [0.09 0.47 0.58; 0.78 0.34 0.15; 0.32 0.37 0.70; 0.48 0.28 0.51]; % Match the four strip orientations.
figures.atlas = figure(Name="Eight electrode sections",Color="w",Position=[60 70 1420 800]); % Create a page for all eight levels.
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
    fill(ax,x,y,[0.97 0.91 0.77],EdgeColor=[0.35 0.29 0.22],LineWidth=1.5); % Show epineurium and its outer edge.
    for j = 1:size(fasc,1) % Draw the perineurium and endoneurium of each fascicle.
        xyOuter = ellipsePoints(fasc(j,:),p.geometry.tPeri_m,theta); % Get the outer perineurial boundary.
        xyInner = ellipsePoints(fasc(j,:),0,theta); % Get the endoneurial boundary.
        fill(ax,mm*xyOuter(:,1),mm*xyOuter(:,2),[0.65 0.40 0.35],EdgeColor="none"); % Show the perineurium as a distinct band.
        fill(ax,mm*xyInner(:,1),mm*xyInner(:,2),[0.96 0.75 0.68],EdgeColor=[0.53 0.29 0.28]); % Show endoneurium within that band.
        text(ax,mm*fasc(j,1),mm*fasc(j,2),string(j),HorizontalAlignment="center",FontSize=8); % Identify each fascicle by its current row number.
    end % Finish the ten provisional fascicles.
    path = mm*strip.path_m(:,1:2); % Convert the routed strip centreline to millimetres.
    plot(ax,path(:,1),path(:,2),Color=colors(abs(strip.section),:),LineWidth=2); % Show the route that avoids fascicles.
    for j = 1:size(strip.padFaces_m,3) % Show all four internal active faces at this level.
        face = mm*strip.padFaces_m(:,:,j); % Read the finite pad vertices.
        plot(ax,face(1:2,1),face(1:2,2),Color=[0.06 0.14 0.20],LineWidth=5); % Project the active pad width.
    end % Finish the internal pads.
    cuff = model.nerve.cuffFaces_m; % Read all external pad faces.
    for j = 1:size(cuff,3) % Find the two cuff contacts at this station.
        face = mm*cuff(:,:,j); % Convert cuff vertices to millimetres.
        if abs(mean(face(:,3))-mm*z_m)<1e-5 % Select the current axial level.
            plot(ax,face(1:2,1),face(1:2,2),Color=[0.05 0.11 0.15],LineWidth=6); % Show the opposing cuff pad.
        end % Skip contacts at other levels.
    end % Finish locating cuff contacts.
    axis(ax,'equal'); % Preserve anatomical proportions.
    xlim(ax,mm*max(p.geometry.nerveSemiAxes_m)*[-1.12 1.12]); % Use the same x scale in all panels.
    ylim(ax,mm*max(p.geometry.nerveSemiAxes_m)*[-1.12 1.12]); % Use the same y scale in all panels.
    title(ax,sprintf('Section %+d  |  z = %.1f mm',strip.section,mm*z_m)); % Give the physical location and station ID.
    xlabel(ax,'x [mm]'); % Label transverse x position.
    ylabel(ax,'y [mm]'); % Label transverse y position.
end % Finish the atlas of all stations.
title(layout,'All eight stations | endoneurium pink, perineurium brown, epineurium gold'); % Explain the tissue colours.
%% 3. Highlight one real detour around the central fascicle
figures.deviation = figure(Name="Strip deflection detail",Color="w",Position=[130 120 870 760]); % Create a large explanation view.
ax = axes(figures.deviation); % Create axes for a single section.
hold(ax,'on'); % Overlay anatomy, nominal path and routed path.
strip = model.strips(1); % Select section -4 for the close-up.
[~,iz] = min(abs(model.nerve.z_m(:,1)-strip.path_m(1,3))); % Select its flattened cross-section.
fill(ax,mm*model.nerve.x_m(iz,:),mm*model.nerve.y_m(iz,:),[0.97 0.91 0.77], ...
    EdgeColor=[0.35 0.29 0.22],LineWidth=1.6); % Show the epineurium and outer edge.
for j = 1:size(fasc,1) % Draw the separate tissue sheaths for every fascicle.
    outer = ellipsePoints(fasc(j,:),p.geometry.tPeri_m,theta); % Compute the perineurial exterior.
    inner = ellipsePoints(fasc(j,:),0,theta); % Compute the endoneurial interior.
    fill(ax,mm*outer(:,1),mm*outer(:,2),[0.65 0.40 0.35],EdgeColor="none"); % Display the perineurial band.
    fill(ax,mm*inner(:,1),mm*inner(:,2),[0.96 0.75 0.68],EdgeColor=[0.53 0.29 0.28]); % Display the endoneurium.
    text(ax,mm*fasc(j,1),mm*fasc(j,2),"F"+j,HorizontalAlignment="center"); % Label fascicles without assumed functional assignments.
end % Finish anatomical layers.
path = mm*strip.path_m(:,1:2); % Read the accepted route in millimetres.
routedLine = plot(ax,path(:,1),path(:,2),Color=[0.05 0.40 0.62],LineWidth=3,DisplayName="Routed strip"); % Trace the route outside the perineurium.
straightLine = plot(ax,[path(1,1),path(end,1)],[path(1,2),path(end,2)],'--', ...
    Color=[0.78 0.19 0.17],LineWidth=1.8,DisplayName="Straight reference"); % Compare with the nominal straight route.
for j = 1:size(model.nerve.cuffFaces_m,3) % Show the external contacts at this selected section.
    cuffFace = mm*model.nerve.cuffFaces_m(:,:,j); % Convert this cuff face to millimetres.
    if abs(mean(cuffFace(:,3))-mm*strip.path_m(1,3))<1e-5 % Select this section's opposing cuff pair.
        plot(ax,cuffFace(1:2,1),cuffFace(1:2,2),Color=[0.02 0.10 0.16],LineWidth=6); % Mark the rigid contact surface.
    end % Skip the other axial stations.
end % Finish displaying the external contacts.
for j = 1:size(strip.padFaces_m,3) % Show each internal contact on its routed position.
    face = mm*strip.padFaces_m(:,:,j); % Read this finite active face.
    plot(ax,face(1:2,1),face(1:2,2),Color=[0.02 0.10 0.16],LineWidth=6); % Highlight the contact width.
end % Finish drawing internal active sites.
axis(ax,'equal'); % Keep the comparison anatomically proportionate.
lim = mm*max(p.geometry.nerveSemiAxes_m)*1.12; % Frame the full nerve section.
xlim(ax,[-lim lim]); ylim(ax,[-lim lim]); % Use equal transverse viewing limits.
xlabel(ax,'x [mm]'); ylabel(ax,'y [mm]'); % State the position units.
title(ax,sprintf('Section %+d | geometrical detour %.3f mm',strip.section,mm*abs(strip.offset_m))); % Report the planned displacement without claiming tissue mechanics.
legend(ax,[routedLine,straightLine],{'Routed strip','Straight reference'},'Location','southoutside'); % Explain only the two compared routes.
end % Return both anatomy figures.
function xy = ellipsePoints(f,peri_m,theta) % Sample one rotated fascicle and optional sheath.
a = f(3)+peri_m; b = f(4)+peri_m; % Read both semi-axes with the outer perineurial thickness.
c = cosd(f(5)); s = sind(f(5)); % Rotate the fascicle in the transverse plane.
xy = [f(1)+a*cos(theta(:))*c-b*sin(theta(:))*s, ...
      f(2)+a*cos(theta(:))*s+b*sin(theta(:))*c]; % Return boundary coordinates in metres.
end % Finish the ellipse sampler.

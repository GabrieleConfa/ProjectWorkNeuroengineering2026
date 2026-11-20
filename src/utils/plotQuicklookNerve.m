function fig = plotQuicklookNerve(model) % Draw the same device geometry in 3D and along the nerve axis.
%% 1. Prepare the MATLAB views
mm = 1e3; % Convert stored coordinates from metres to millimetres.
p = model.parameters; % Read the shared nerve and drawing inputs.
nerve = model.nerve; % Read the sampled nerve surface and cuff polygons.
fasc = p.geometry.fascicles_m; % Read the straight fascicle axes and radii.
fig = figure(Name="Median nerve 3D | " + model.mode, Color="w", Position=[60 110 1190 670]); % Create the interactive figure.
layout = tiledlayout(fig,1,2,TileSpacing="compact",Padding="compact"); % Show both views side by side.
space = nexttile(layout); % Allocate the oblique 3D view.
hold(space,'on'); % Overlay nerve, fascicles and contact faces.
front = nexttile(layout); % Allocate the end-on projection.
hold(front,'on'); % Overlay the transverse stations at their projected locations.
colors = [0.12 0.49 0.60; 0.77 0.35 0.16; 0.34 0.38 0.69; 0.50 0.30 0.51]; % Identify the four mirrored station pairs.
theta = linspace(0,2*pi,181); % Sample the fascicle cross-sections.
z_m = p.domain.length_m/2 * [-1; 1]; % Use the full axial extent for each fascicle.

%% 2. Draw the nerve and fascicles
surf(space, mm*nerve.x_m, mm*nerve.y_m, mm*nerve.z_m, ...
    FaceColor=[0.76 0.65 0.49], FaceAlpha=0.14, EdgeColor="none"); % Draw the locally flattened outer surface.
plot(front, mm*nerve.x_m(1,:), mm*nerve.y_m(1,:), ...
    Color=[0.26 0.27 0.32], LineWidth=1.5); % Plot the nerve silhouette beyond the electrode stations.
for k = 1:size(fasc,1) % Draw the fixed outer perineurium boundaries.
    a_m = fasc(k,3) + p.geometry.tPeri_m; % Include perineurium along the first fascicle axis.
    b_m = fasc(k,4) + p.geometry.tPeri_m; % Include perineurium along the second fascicle axis.
    angle_deg = fasc(k,5); % Read the fascicle rotation around the nerve axis.
    x_mm = mm*(fasc(k,1) + a_m*cos(theta)*cosd(angle_deg) - b_m*sin(theta)*sind(angle_deg)); % Rotate the oval x samples.
    y_mm = mm*(fasc(k,2) + a_m*cos(theta)*sind(angle_deg) + b_m*sin(theta)*cosd(angle_deg)); % Rotate the oval y samples.
    surf(space, repmat(x_mm,2,1), repmat(y_mm,2,1), mm*repmat(z_m,1,numel(theta)), ...
        FaceColor=[0.70 0.43 0.31], FaceAlpha=0.26, EdgeColor="none"); % Draw the fascicle along the nerve axis.
    patch(front, x_mm, y_mm, [0.97 0.88 0.80], EdgeColor=[0.60 0.44 0.37]); % Show the same fascicle in transverse projection.
    text(front, mm*fasc(k,1), mm*fasc(k,2), string(k), ...
        HorizontalAlignment="center", Color=[0.31 0.23 0.22], FontSize=8); % Number the fascicles for overlap review.
end % Finish drawing the anatomy.

%% 3. Draw the pad faces and transverse stations
for k = 1:numel(model.strips) % Draw each station at its measured z coordinate.
    strip = model.strips(k); % Read the route and finite pad geometry.
    c = colors(abs(strip.section),:); % Match mirrored sections by colour.
    path = mm*strip.path_m; % Convert the strip centreline for plotting.
    plot3(space,path(:,1),path(:,2),path(:,3),Color=c,LineWidth=1.4); % Show the strip path through the nerve.
    if strip.section < 0 % Draw each coincident mirrored pair only once in the end-on projection.
        plot(front,path(:,1),path(:,2),Color=c,LineWidth=1.2); % Overlay that pair's transverse route.
    end % Finish drawing the projection of the strip.
    for j = 1:size(strip.padFaces_m,3) % Draw each internal active site separately.
        drawPad(space, front, mm*strip.padFaces_m(:,:,j), c, strip.section < 0); % Show the ideal contact patch and its projected interval.
    end % Finish plotting this station's internal pads.
end % Finish plotting the internal part of the device.
for k = 1:size(nerve.cuffFaces_m,3) % Draw the outer contacts on the flattened nerve.
    face = mm*nerve.cuffFaces_m(:,:,k); % Convert this finite active face to millimetres.
    drawPad(space, front, face, [0.14 0.22 0.26], mean(face(:,3)) < 0); % Draw the cuff face in the nerve's flat plane.
end % Finish plotting the external contacts.

%% 4. Set equal transverse scales and label the two views
axis(space,'equal'); % Keep the 3D geometry in the same units along every axis.
view(space, 44, 26); % Start with an oblique view of the nerve and all stations.
grid(space,'on'); % Show the axial spacing against the millimetre scale.
xlabel(space,'x [mm]'); % Label the transverse x coordinate.
ylabel(space,'y [mm]'); % Label the transverse y coordinate.
zlabel(space,'z [mm]'); % Label the nerve axis and station positions.
title(space,'3D | 8 electrode stations'); % Identify the full device view.
axis(front,'equal'); % Keep the front view circular in the undeformed region.
xlim(front, mm*max(p.geometry.nerveSemiAxes_m)*[-1.1 1.1]); % Frame the transverse geometry.
ylim(front, mm*max(p.geometry.nerveSemiAxes_m)*[-1.1 1.1]); % Use the same view extent in both directions.
xlabel(front,'x [mm]'); % Label the frontal x coordinate.
ylabel(front,'y [mm]'); % Label the frontal y coordinate.
title(front, ["End-on projection | 4 mirrored pairs", "Thin lines: strips  •  Thick lines: active faces"]); % Explain the overlaid station pairs and contacts.
title(layout,"Median nerve | " + model.mode); % Label the interaction scenario.
end % Return the interactive figure.

function drawPad(space, front, vertices, color, showProjection) % Draw a single finite active surface.
patch(space,vertices(:,1),vertices(:,2),vertices(:,3),color,EdgeColor="none"); % Render the electrode face without a bulk electrode volume.
if showProjection % Project one face per mirrored pair to the transverse plane.
    plot(front,vertices(1:2,1),vertices(1:2,2),Color=color,LineWidth=4); % Make the projected active interval visible.
end % Finish drawing the projected site.
end % Finish drawing the pad.

function fig = plotFascicleMap(p) % Plot the provisional anatomy and the group's editable functional assignments.
%% 1. Anatomical map at a station outside the cuff
f = p.geometry.fascicles_m; % Read centres, semiaxes and rotations from the current case.
labels = repmat("pending",size(f,1),1); % Keep missing functions explicit.
if isfield(p.geometry,'functionLabels'), labels=string(p.geometry.functionLabels(:)); end % Accept one assignment per fascicle.
assert(numel(labels)==size(f,1),'Supply one function label per fascicle.'); % Reject a mismatched mapping.
fig = figure(Name="Fascicle dimensions and functions",Color="w",Theme="light",Position=[60 60 1300 600]); % Create a readable anatomical figure.
layout = tiledlayout(fig,1,2); % Place the measurement table next to the map.
ax = nexttile(layout); % Select the anatomical panel.
angles = linspace(0,2*pi,361)'; % Close the nerve and fascicle outlines.
ray = [cos(angles),sin(angles)]; % Define the sampled planar directions.
xy = 1e3*ray.*nerveBoundaryRadius(p,ray); % Convert the uncompressed outline to millimetres.
fill(ax,xy(:,1),xy(:,2),[0.97 0.91 0.77]); % Colour the epineurium.
hold(ax,'on'); % Add each fascicle to the same anatomical section.
for k = 1:size(f,1) % Draw and identify each provisional structure.
    rot = [cosd(f(k,5)) sind(f(k,5));-sind(f(k,5)) cosd(f(k,5))]; % Rotate the ellipse into anatomical coordinates.
    for thickness = [p.geometry.tPeri_m,0] % Draw the sheath before its endoneurial interior.
        edge = 1e3*(f(k,1:2)+(ray.*(f(k,3:4)+thickness))*rot); % Compute the two tissue outlines.
        colour = [0.96 0.75 0.68]; % Set the endoneurium colour.
        if thickness>0, colour=[0.65 0.40 0.35]; end % Use the same sheath colour as the station atlas.
        fill(ax,edge(:,1),edge(:,2),colour,EdgeColor='none'); % Draw this tissue layer.
    end % Finish the fascicle outline.
    text(ax,1e3*f(k,1),1e3*f(k,2),"F"+k,HorizontalAlignment='center'); % Match the labels in the solver.
end % Finish drawing all fascicles.
axis(ax,'equal'); xlabel(ax,'x (mm)'); ylabel(ax,'y (mm)'); % Keep the section to scale.
title(ax,sprintf('Perineurium %.1f µm | provisional anatomy',1e6*p.geometry.tPeri_m)); % Display the configured thickness.
%% 2. List diameters and functional assignments
ax = nexttile(layout); axis(ax,'off'); % Reserve the right panel for the numerical key.
headers=["ID","x (mm)","y (mm)","D1 (mm)","D2 (mm)","Angle (°)","Function"]; % Distinguish full diameters from the stored semiaxes.
data=["F"+(1:size(f,1))',compose("%.3f",1e3*f(:,1:2)),compose("%.3f",2e3*f(:,3:4)),compose("%.0f",f(:,5)),labels]; % Populate the table from the same parameters as the map.
cols=[0.05 0.15 0.29 0.43 0.57 0.71 0.86]; % Align the identifiers, geometry and functional labels by column.
nRows=size(data,1); axis(ax,[0 1 0 nRows+1]); % Give each fascicle its own row.
hold(ax,'on'); % Compose a table that exports with standard MATLAB graphics.
for row=0:nRows % Draw the header and the provisional data rows.
    y=nRows-row; colour=[0.96 0.97 0.98]; txt=[0.12 0.17 0.23]; values=headers; % Start with the header labels.
    if row==0, colour=[0.15 0.24 0.32]; txt=[1 1 1]; else, values=data(row,:); if mod(row,2)==0, colour=[1 1 1]; end; end % Shade alternate fascicles for readability.
    patch(ax,[0 1 1 0],[y y y+1 y+1],colour,EdgeColor='none'); % Fill one full table row.
    for col=1:numel(cols) % Print each value under its corresponding header.
        text(ax,cols(col),y+0.5,values(col),Color=txt,HorizontalAlignment="center",FontSize=9,Interpreter="none"); % Keep numeric precision and function names explicit.
    end % Finish this row.
end % Finish the parameter-driven table.
title(ax,'Provisional geometry and functional assignments'); % Signal that dimensions and positions await group measurements.
end % Return the map and its matching table.

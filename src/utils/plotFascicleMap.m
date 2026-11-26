function fig = plotFascicleMap(p) % Plot the anatomy of one section and the table of fascicle dimensions.
%% 1. Anatomical map at a station outside the cuff
f = p.geometry.fascicles_m; % Read centres, semiaxes and rotations from the current case.
fig = figure(Name="Fascicle map and dimensions",Color="w",Theme="light",Position=[60 60 1300 600]); % Create a readable anatomical figure.
layout = tiledlayout(fig,1,2); % Place the measurement table next to the map.
ax = nexttile(layout); % Select the anatomical panel.
angles = linspace(0,2*pi,361)'; % Close the nerve and fascicle outlines.
ray = [cos(angles),sin(angles)]; % Define the sampled planar directions.
xy = 1e3*ray.*nerveBoundaryRadius(p,ray); % Convert the uncompressed outline to millimetres.
fill(ax,xy(:,1),xy(:,2),[0.97 0.91 0.77]); % Colour the epineurium.
hold(ax,'on'); % Add each fascicle to the same anatomical section.
polys = stationFascicles(p,0); % Undisturbed outlines (before insertion).
drawFascicles(ax,polys,p.geometry.tPeri_m,1e3,true); % Draw epineurium-free fascicle outlines with perineurium.
axis(ax,'equal'); xlabel(ax,'x (mm)'); ylabel(ax,'y (mm)'); % Keep the section to scale.
title(ax,sprintf('Perineurium %.1f µm',1e6*p.geometry.tPeri_m)); % Display the configured thickness.
%% 2. Table of fascicle dimensions
ax = nexttile(layout); axis(ax,'off'); % Reserve the right panel for the numerical key.
headers=["ID","x (mm)","y (mm)","Area (mm²)","Round.","Feret (mm)"]; % Shape descriptors valid for irregular outlines.
stats=zeros(numel(polys),5); for k=1:numel(polys), P=polys{k}; feret=max(max(sqrt((P(:,1)-P(:,1)').^2+(P(:,2)-P(:,2)').^2))); A=polyarea(P(:,1),P(:,2)); stats(k,:)=[mean(P,1)*1e3,A*1e6,4*A/(pi*feret^2),feret*1e3]; end % Centroid, area, roundness 4A/(pi*Feret^2) and maximum Feret diameter.
data=["F"+(1:numel(polys))',compose("%.3f",stats(:,1:2)),compose("%.3f",stats(:,3)),compose("%.2f",stats(:,4)),compose("%.2f",stats(:,5))]; % Populate the table from the same outlines as the map.
cols=[0.08 0.24 0.40 0.57 0.74 0.90]; % Six centred columns, wide enough not to overlap.
nRows=size(data,1); axis(ax,[0 1 0 nRows+1]); % Give each fascicle its own row.
hold(ax,'on'); % Compose a table that exports with standard MATLAB graphics.
for row=0:nRows % Draw the header and the provisional data rows.
    y=nRows-row; colour=[0.96 0.97 0.98]; txt=[0.12 0.17 0.23]; values=headers; % Start with the header labels.
    if row==0, colour=[0.15 0.24 0.32]; txt=[1 1 1]; else, values=data(row,:); if mod(row,2)==0, colour=[1 1 1]; end; end % Shade alternate fascicles for readability.
    patch(ax,[0 1 1 0],[y y y+1 y+1],colour,EdgeColor='none'); % Fill one full table row.
    for col=1:numel(cols) % Print each value under its corresponding header.
        text(ax,cols(col),y+0.5,values(col),Color=txt,HorizontalAlignment="center",FontSize=8.5,Interpreter="none"); % Keep numeric precision explicit.
    end % Finish this row.
end % Finish the parameter-driven table.
title(ax,'Fascicle dimensions');
end % Return the map and its matching table.

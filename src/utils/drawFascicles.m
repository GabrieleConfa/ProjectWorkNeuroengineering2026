function h = drawFascicles(ax, polys, tPeri, scale, labels) % Draw perineurium (band) and endoneurium of fascicle outlines.
%% Purpose
% Shared drawing of fascicles for all figures, valid for circles, ellipses and random irregular outlines.
% polys: cell of Nx2 outlines (metres). scale: unit factor (1e3 for millimetres). labels: draw "F<k>" at the centroid.
if nargin<4, scale = 1e3; end, if nargin<5, labels = true; end
hold(ax,'on'); h = struct('peri',[],'endo',[]); % Handles for legends.
for k = 1:numel(polys)
    P = polys{k}; outer = polybuffer(polyshape(P),tPeri); % Outer perineurial boundary.
    h.peri = fill(ax,scale*outer.Vertices(:,1),scale*outer.Vertices(:,2),[0.65 0.40 0.35],EdgeColor="none"); % Perineurium band.
    h.endo = fill(ax,scale*P(:,1),scale*P(:,2),[0.96 0.75 0.68],EdgeColor=[0.53 0.29 0.28]); % Endoneurium.
    if labels, text(ax,scale*mean(P(:,1)),scale*mean(P(:,2)),"F"+k,HorizontalAlignment="center",FontSize=8); end % Identify the fascicle.
end
end

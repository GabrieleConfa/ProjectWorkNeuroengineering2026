function drawPersonalStation(g) % Draw geometry from the actual 2D solution in millimetres.
outline=g.outline([1:end 1],:)*1e3; plot(outline(:,1),outline(:,2),'k-',LineWidth=0.9); % Follow both rigid contact footprints.
th=linspace(0,2*pi,241)'; circle=[cos(th) sin(th)]; % Sample closed tissue boundaries.
for k=1:numel(g.polys) % Outer perineurial boundary of each fascicle (circles, ellipses or random irregular shapes).
    xy=1e3*polybuffer(polyshape(g.polys{k}),g.q.geometry.tPeri_m).Vertices; xy=[xy;xy(1,:)]; % Close the outline.
    plot(xy(:,1),xy(:,2),'k-',LineWidth=0.6); % Mark the tissue boundary on each field.
end % Finish the sheaths.
plot(1e3*g.path(:,1),1e3*g.path(:,2),Color=[0.13 0.38 0.92],LineWidth=1.1); % Show the continuous carrier as a geometry overlay.
end

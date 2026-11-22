function drawPersonalStation(g) % Draw geometry from the actual 2D solution in millimetres.
outline=g.outline([1:end 1],:)*1e3; plot(outline(:,1),outline(:,2),'k-',LineWidth=0.9); % Follow both rigid contact footprints.
th=linspace(0,2*pi,241)'; circle=[cos(th) sin(th)]; % Sample closed tissue boundaries.
for k=1:size(g.f,1) % Use both semiaxes and rotation for each outer perineurium.
    f=g.f(k,:); rotate=[cosd(f(5)) sind(f(5));-sind(f(5)) cosd(f(5))]; % Rotate into anatomical coordinates.
    xy=1e3*(f(1:2)+(circle.*(f(3:4)+g.q.geometry.tPeri_m))*rotate); % Project the configured outer sheath.
    plot(xy(:,1),xy(:,2),'k-',LineWidth=0.6); % Mark the tissue boundary on each field.
end % Finish the ten sheaths.
plot(1e3*g.path(:,1),1e3*g.path(:,2),Color=[0.13 0.38 0.92],LineWidth=1.1); % Show the continuous carrier as a geometry overlay.
end

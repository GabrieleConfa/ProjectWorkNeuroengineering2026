function radius_m = nerveBoundaryRadius(p, directions) % Find an elliptical nerve boundary along unit rays.
axes_m = p.geometry.nerveSemiAxes_m; % Read the editable transverse semi-axes.
angle_deg = p.geometry.nerveRotation_deg; % Read the editable rotation of the ellipse.
assert(numel(axes_m)==2 && all(axes_m>0), 'Nerve semi-axes must be positive.'); % Reject missing anatomical sizes.
assert(size(directions,2)==2, 'Directions must have x and y columns.'); % Require transverse unit directions.
c = cosd(angle_deg); % Compute the major-axis rotation cosine.
s = sind(angle_deg); % Compute the major-axis rotation sine.
localX = c*directions(:,1) + s*directions(:,2); % Rotate ray directions into the nerve frame.
localY = -s*directions(:,1) + c*directions(:,2); % Compute their local minor-axis component.
radius_m = 1 ./ sqrt((localX/axes_m(1)).^2 + (localY/axes_m(2)).^2); % Intersect each ray with the ellipse.
end % Return radial boundary distances in metres.

function [x,y] = fascicleOutline(F,k,theta) % Sample a supplied ellipse or a legacy circular fascicle in metres.
a = F.rEndo_m(k); b = a; phi = 0; % Read the circular reference anatomy.
if ismember('a_m',F.Properties.VariableNames) % Detect editable elliptical anatomy.
    a = F.a_m(k); b = F.b_m(k); phi = F.angle_deg(k); % Read both axes and their orientation.
end % Finish resolving the input representation.
x = F.x_m(k)+a*cos(theta)*cosd(phi)-b*sin(theta)*sind(phi); % Rotate and translate the x coordinate.
y = F.y_m(k)+a*cos(theta)*sind(phi)+b*sin(theta)*cosd(phi); % Rotate and translate the y coordinate.
end % Return the physical outline in the same frame as the mesh.

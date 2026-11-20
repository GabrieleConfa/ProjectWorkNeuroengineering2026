function meshPreview = buildQuicklookMesh(p, cuffSites) % Mesh a short 3D slab through one flattened electrode station.
%% 1. Define the station and its 2D material regions
assert(all(cuffSites.section == cuffSites.section(1)), 'Choose one cuff pair from one station.'); % Keep a single axial slab.
assert(all(abs(p.geometry.nerveSemiAxes_m - p.geometry.radius_m) < 1e-12) && ...
    all(abs(p.geometry.fascicles_m(:,3)-p.geometry.fascicles_m(:,4)) < 1e-12), ...
    'The local slab mesher currently requires circular nerve and fascicles.'); % Avoid silently meshing the wrong oval anatomy.
R = p.geometry.radius_m; % Read the provisional nerve radius.
depth = p.electrode.cuffIndentation_m; % Read the prescribed surface flattening.
angles = linspace(0, 2*pi, 241); % Sample the undeformed circular nerve boundary.
xy = R*[cos(angles(:)), sin(angles(:))]; % Start with a circular cross-section.
normals = [cosd(cuffSites.theta_deg), sind(cuffSites.theta_deg)]; % Read the two rigid pad normals.
for k = 1:height(cuffSites) % Clip the nerve to each rigid cuff plane.
    projection = xy*normals(k,:)'; % Find the signed distance along this pad normal.
    cut = projection > R-depth; % Find sampled points beyond the cuff plane.
    xy(cut,:) = xy(cut,:) .* ((R-depth) ./ projection(cut)); % Pull those boundary points onto the flat plane.
end % Finish clipping the nerve section.
xy(end,:) = []; % Remove the repeated first vertex.
numFasc = size(p.geometry.fascicles_m,1); % Count the provisional fascicles.
numShapes = 1 + 2*numFasc; % Reserve one nerve polygon and two circles per fascicle.
maxRows = max(2*size(xy,1)+2, 10); % Size the geometry matrix for the boundary polygon.
gd = zeros(maxRows,numShapes); % Allocate the MATLAB constructive geometry matrix.
gd(1:2*size(xy,1)+2,1) = [2; size(xy,1); xy(:,1); xy(:,2)]; % Encode the flattened nerve boundary.
names = strings(1,numShapes); % Allocate the region names.
names(1) = "N"; % Name the outer nerve region.
for k = 1:numFasc % Add endoneurium and perineurium circles.
    f = p.geometry.fascicles_m(k,:); % Read one fascicle centre and endoneurium radius.
    gd(1:4,2*k) = [1; f(1); f(2); f(3)+p.geometry.tPeri_m]; % Encode its outer perineurial boundary.
    gd(1:4,2*k+1) = [1; f(1); f(2); f(3)]; % Encode its endoneurium boundary.
    names(2*k:2*k+1) = ["P"+k, "E"+k]; % Give the two nested regions distinct names.
end % Finish constructing tissue regions.
[dl,bt] = decsg(gd,char(strjoin(names,"+")),char(names(:))'); % Decompose all regions into nonoverlapping faces.
geometry2D = fegeometry(dl); % Create a 2D PDE geometry for extrusion.
%% 2. Extrude the section and create a tetrahedral mesh
height_m = cuffSites.sideZ_m(1); % Match the axial length of the active pad.
geometry3D = extrude(geometry2D,height_m); % Form an axial volume with conforming tissue interfaces.
geometry3D = translate(geometry3D,[0,0,cuffSites.z_m(1)-height_m/2]); % Place the slab at its measured station.
geometry3D = generateMesh(geometry3D,Hmax=0.12e-3,Hmin=0.02e-3,GeometricOrder="linear"); % Create an initial tetrahedral mesh.
mesh = geometry3D.Mesh; % Extract the tetrahedral node and element arrays.
meshPreview.geometry = geometry3D; % Return the labelled volume geometry.
meshPreview.mesh = mesh; % Return its tetrahedra and node coordinates.
meshPreview.regionTable = bt; % Retain the 2D membership matrix for later material assignment.
meshPreview.station = cuffSites.section(1); % Record the station used for this local mesh.
end % Return the local 3D mesh preview.

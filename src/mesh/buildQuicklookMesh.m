function meshPreview = buildQuicklookMesh(p, cuffSites) % Mesh one contact slab using the same local contour as the preview.
%% 1. Define the station and its 2D material regions
assert(all(cuffSites.section == cuffSites.section(1)), 'Choose one cuff pair from one station.'); % Keep a single axial slab.
angles = linspace(0,2*pi,241)'; % Sample the nerve outline around one pad-centred section.
angles(end) = []; % Remove the duplicate vertex at the polygon seam.
directions = [cos(angles),sin(angles)]; % Define rays for either a circular or oval nerve.
radii_m = cuffDeformedRadius(p,cuffSites,cuffSites.z_m(1),directions); % Reuse the local geometric contact model.
xy = radii_m(:).*directions; % Form the actual pad-contacted cross-section.
numFasc = size(p.geometry.fascicles_m,1); % Count the provisional fascicles.
numShapes = 1 + 2*numFasc; % Reserve the nerve polygon and two ellipses per fascicle.
maxRows = max(2*size(xy,1)+2, 10); % Size the geometry matrix for the boundary polygon.
gd = zeros(maxRows,numShapes); % Allocate the MATLAB constructive geometry matrix.
gd(1:2*size(xy,1)+2,1) = [2; size(xy,1); xy(:,1); xy(:,2)]; % Encode the flattened nerve boundary.
names = strings(1,numShapes); % Allocate the region names.
names(1) = "N"; % Name the outer nerve region.
for k = 1:numFasc % Add endoneurium and perineurium ellipses.
    f = p.geometry.fascicles_m(k,:); % Read one fascicle centre and endoneurium radius.
    gd(1:6,2*k) = [4; f(1); f(2); f(3)+p.geometry.tPeri_m; f(4)+p.geometry.tPeri_m; deg2rad(f(5))]; % Encode the rotated outer sheath.
    gd(1:6,2*k+1) = [4; f(1); f(2); f(3); f(4); deg2rad(f(5))]; % Encode the rotated endoneurium.
    names(2*k:2*k+1) = ["P"+k, "E"+k]; % Give the two nested regions distinct names.
end % Finish constructing tissue regions.
[dl,bt] = decsg(gd,char(strjoin(names,"+")),char(names(:))'); % Decompose all regions into nonoverlapping faces.
geometry2D = fegeometry(dl); % Create a 2D PDE geometry for extrusion.
%% 2. Extrude the section and create a tetrahedral mesh
height_m = cuffSites.sideZ_m(1); % Match the axial length of the active pad.
geometry3D = extrude(geometry2D,height_m); % Form an axial volume with conforming tissue interfaces.
geometry3D = translate(geometry3D,[0,0,cuffSites.z_m(1)-height_m/2]); % Place the slab at its measured station.
geometry3D = generateMesh(geometry3D,Hmax=p.meshPreview.hMax_m,Hmin=p.meshPreview.hMin_m,GeometricOrder="linear"); % Apply editable mesh sizes.
mesh = geometry3D.Mesh; % Extract the tetrahedral node and element arrays.
meshPreview.geometry = geometry3D; % Return the labelled volume geometry.
meshPreview.mesh = mesh; % Return its tetrahedra and node coordinates.
meshPreview.regionTable = bt; % Retain the 2D membership matrix for later material assignment.
meshPreview.station = cuffSites.section(1); % Record the station used for this local mesh.
end % Return the local 3D mesh preview.

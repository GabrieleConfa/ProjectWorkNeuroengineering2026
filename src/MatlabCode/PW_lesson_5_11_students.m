%% =====================================================================
%  PERIPHERAL NERVE MODEL - PART 1
%  Geometry, tissue electrical properties and mesh of the nerve cross-section
%
%  Requirements: MATLAB R2023a or later + Partial Differential Equation Toolbox
%  Units:        lengths in mm, conductivities in S/m
%
%  Exercises 0-3 completed: tissue values in section 2 are from the literature
%  (Davis et al. 2023); see the comments next to each value.
%
%  Workflow:
%    1. Geometry parameters          6. Tissue of each geometry face
%    2. Tissue electrical properties 7. Assign material properties
%    3. Mesh parameters              8. Generate the mesh
%    4. Build the 2D geometry        9. Link mesh elements to tissues (+ check)
%    5. Create the FE model         10-14. Plots
%
%  General idea of a finite element (FE) simulation:
%    - the domain (here the nerve cross-section) is split into many small,
%      simple pieces called ELEMENTS (here triangles), connected at NODES;
%    - each element:
%        * represents a small piece of the geometry, defined by its nodes
%        * carries the material properties of the tissue it belongs to
%        * approximates the solution inside it from the values at its nodes
%        * contributes its own equations to the global system to be solved
%    - this script prepares everything up to the mesh: the solution
%      (boundary conditions, sources, solve) is Part 2.
%  =====================================================================

clear; close all; clc;       % clear variables, close figures, clear the Command Window


%% ---------- 1. GEOMETRY PARAMETERS ----------
% All the numbers that define the shape of the model are collected here, so
% that the geometry can be changed without touching the rest of the code.
%
% The model is a 2D CROSS-SECTION of the nerve (a "slice" perpendicular to its
% axis). From the outside to the inside we find:
%   - surrounding medium : the tissue around the nerve (e.g. muscle)
%   - epineurium         : connective tissue that wraps the whole nerve
%   - perineurium        : thin sheath around each fascicle (very resistive)
%   - endoneurium        : inside each fascicle, where the nerve fibers run
%
% No variable needs to be "declared" in MATLAB: a variable is created the
% first time we assign a value to it.


% ======================= EXERCISE 0  =======================

Rnerve  = 1.5;            % nerve radius [mm]
                          % = radius of a circle with the area of a proximal-forearm
                          %   median nerve, 7.13 mm^2 (Delgado-Martinez 2016, p. 6).
                          %   Must be > ~1.29 mm, or the outer fascicles stick out.
Rmedium = 3 * Rnerve;     % radius of the surrounding medium [mm]
                          % (large enough that its outer edge is "far" from the nerve)
tPeri   = 0.04;           % perineurium thickness, same for all fascicles [mm]
                          % = 5.5-8% of these fascicle diameters: on the thick side
                          %   (Grinberg 2008: 3%; Davis 2023 human-vagus fit: 28-37 um).
                          %   If changed, change hPeri too (>= 2 elements across).


%% FASCICLES SETUP
% Fascicles: one row per fascicle  ->  [x_center  y_center  radius]  [mm]
% "radius" is the endoneurium radius; the perineurium is added around it.
% The geometry is FIXED (not random): the positions were chosen so that
% fascicles do not overlap each other or the nerve boundary.
fasc = [ 0.00   0.00   0.36
         0.86   0.40   0.30
        -0.84   0.42   0.30
         0.10  -0.92   0.32
        -0.74  -0.58   0.24
         0.86  -0.52   0.24
         0.02   0.96   0.26 ];
nFasc = size(fasc, 1);    % size(A,1) = number of rows = number of fascicles

Lnerve = 8;               % length of the nerve segment drawn in the 3D view [mm]
                          % (used only for the picture, not part of the model)
                          %5 to 8 cm nerve lenght in the simulation, big
                          %device.


%% ---------- 2. TISSUE ELECTRICAL PROPERTIES  ----------
%  Every tissue is described by a few electrical properties. They are stored in
% a TABLE: one row per tissue, one column per property. A table is easy to
% read, to print and to access by name (e.g. tissues.sigmaT).
%
% Tissue numbering used in the whole script:
%   1 = endoneurium, 2 = perineurium, 3 = epineurium, 4 = surrounding medium
%
% Properties:
%   sigmaT = transverse conductivity, i.e. in the cross-section plane   [S/m]
%   sigmaL = longitudinal conductivity, i.e. along the nerve axis       [S/m]
%   epsR   = relative permittivity (used only for the impedance plot)   [-]
%
% Conductivity sigma says how easily current flows through a material
% (high sigma = good conductor). Some tissues are ANISOTROPIC: the fibers make
% current flow more easily along the nerve than across it (sigmaL > sigmaT).
%
% The values below are taken from the literature (source next to them).
 
tissueName = ["Endoneurium"; "Perineurium"; "Epineurium"; "Surrounding medium"];



% ======================= EXERCISE 1 (part a) =======================
% Search the literature for the electrical conductivity of each tissue and
% replace every NaN with a value in S/m.
% Order of the values (same as tissueName):
%     [endoneurium; perineurium; epineurium; surrounding medium]
% Hints:
%   - sigmaT is the conductivity ACROSS the fibers, sigmaL ALONG the fibers
%   - for an ISOTROPIC tissue (same in all directions): sigmaL = sigmaT
%   - the perineurium is much less conductive than the other tissues
%   - for the surrounding medium, assume muscle with fibers parallel to the nerve
%   - write down the source (paper) of every value you use!
% Source: Davis et al. 2023, J Neural Eng 20(2), doi:10.1088/1741-2552/acc42b,
% Table 1 (ASCENT defaults).
% Nerve tissues are given there as resistivities, so sigma = 1/rho:
%   endoneurium rho_T = 6 Ohm*m, rho_L = 1.75 Ohm*m; perineurium 1149 Ohm*m;
%   epineurium 6.3 Ohm*m. Muscle (from Gielen et al. 1984): 0.086 / 0.35 S/m.
sigmaT     = [1/6;     1/1149;  1/6.3;  0.086];   % transverse conductivities   [S/m]
sigmaL     = [1/1.75;  1/1149;  1/6.3;  0.35];    % longitudinal conductivities [S/m]
 
% Automatic check: stops here if some values are still missing
assert(~any(isnan([sigmaT; sigmaL])), ...
       'EXERCISE 1a: replace all the NaN values in sigmaT and sigmaL.');
% ====================================================================
 
tissues = table(tissueName, sigmaT, sigmaL);
 
% Derived quantities, added as new columns of the table.
% NOTE on "./" : the dot means "element by element". 1 ./ v divides 1 by EACH
% element of v and returns a vector of the same size. Without the dot,
% 1 / v would be a MATRIX division (it solves a linear system), which is
% NOT what we want here. Same rule for .* and .^
tissues.rhoT       = 1 ./ tissues.sigmaT;            % transverse resistivity   [Ohm*m]
 
% ======================= EXERCISE 1 (part b) =======================
% Following the example of rhoT (line above), compute two new columns:
%   - rhoL       = longitudinal resistivity [Ohm*m]
%   - anisotropy = ratio between longitudinal and transverse conductivity
%                  (> 1 means that current flows more easily along the nerve)
% Remember to use the element-by-element operators (./) !
tissues.rhoL       = 1 ./ tissues.sigmaL;                % longitudinal resistivity [Ohm*m]
tissues.anisotropy = tissues.sigmaL ./ tissues.sigmaT;   % sigmaL/sigmaT [-] (1 = isotropic)
 
% Automatic check
assert(~any(isnan(tissues.rhoL)) && ~any(isnan(tissues.anisotropy)), ...
       'EXERCISE 1b: compute tissues.rhoL and tissues.anisotropy.');
% ====================================================================
 
disp('Tissue electrical properties (Davis et al. 2023, Table 1):')
disp(tissues)
 
% Colors used for each tissue in the plots, as [red green blue], each 0...1
tissueColor = [0.95 0.70 0.75     % endoneurium
               0.45 0.25 0.35     % perineurium
               0.98 0.90 0.75     % epineurium
               0.85 0.93 1.00];   % surrounding medium
 
 
%% ---------- 3. MESH PARAMETERS ----------
% The mesh splits the geometry into small triangles. Here we choose the target
% size of the triangles (length of their edges) in each tissue [mm].
%
% Rule of thumb:
%   - SMALL elements where the geometry is thin or where we need accuracy
%     (the perineurium is only 0.04 mm thick: we want at least 2 elements
%      across it, otherwise it would not be described correctly);
%   - LARGE elements where we do not need detail (far from the nerve),
%     to keep the number of elements, and so the computation time, low.
 
hEndo   = 0.05;     % endoneurium
hPeri   = 0.02;     % perineurium: half its thickness -> at least 2 elements across
hEpi    = 0.12;     % epineurium
hMedium = 0.60;     % surrounding medium (far from the nerve: large elements)
hGrad   = 1.3;      % growth rate: max size ratio between neighboring elements.
                    % Values close to 1 give smooth (but many) elements.
 
 
%% ---------- 4. 2D GEOMETRY OF THE CROSS-SECTION (decsg) ----------
% Idea: every boundary between two tissues is a circle. We
%   (a) describe all the circles             -> matrix gd
%   (b) give each circle a name              -> ns
%   (c) say how to combine them              -> set formula sf
%   (d) let decsg split the plane into regions (the "faces" of the geometry)
%       -> dl (the geometry) and bt (which circles contain each face)
 
% --- How many circles? ---
% 1 for the surrounding medium (M), 1 for the nerve boundary (N),
% 2 for each fascicle: outer edge of the perineurium (P) and endoneurium (E).
% With 7 fascicles: 2 + 2*7 = 16 circles.
nShapes = 2 + 2*nFasc;
 
% --- Empty "containers", filled below (preallocation) ---
% gd = Geometry Description matrix: ONE COLUMN PER CIRCLE, 4 rows:
%        [ 1        ]  <- shape type: 1 = circle (2 = polygon, 3 = rectangle)
%        [ x_center ]
%        [ y_center ]
%        [ radius   ]
% names = name of each circle, in the same order as the columns of gd.
gd      = zeros(4, nShapes);       % 4 x 16 matrix of zeros
names   = strings(1, nShapes);     % 16 empty strings
 
% --- First two circles: surrounding medium and nerve, both centered in (0,0) ---
% gd(:,1) means "all rows of column 1"
gd(:,1) = [1; 0; 0; Rmedium];   names(1) = "M";
gd(:,2) = [1; 0; 0; Rnerve];    names(2) = "N";
 
% --- Two concentric circles for each fascicle ---
%   P = outer edge of the perineurium: radius = endoneurium radius + tPeri
%   E = endoneurium:                   radius = endoneurium radius
% The indices 1+2*i and 2+2*i place them right after M and N:
%   fascicle 1 -> columns 3 (P1) and 4 (E1)
%   fascicle 2 -> columns 5 (P2) and 6 (E2)   ...   fascicle 7 -> 15 and 16
% "P" + i joins text and number: with i = 2 it gives "P2".
for i = 1:nFasc
    gd(:, 1+2*i) = [1; fasc(i,1); fasc(i,2); fasc(i,3) + tPeri];   names(1+2*i) = "P" + i;
    gd(:, 2+2*i) = [1; fasc(i,1); fasc(i,2); fasc(i,3)];           names(2+2*i) = "E" + i;
end
% Result (first columns):
%            M      N      P1     E1     P2     E2    ...
%   type     1      1      1      1      1      1
%   x        0      0      0      0      0.86   0.86
%   y        0      0      0      0      0.40   0.40
%   r        4.5    1.5    0.40   0.36   0.34   0.30
 
% --- Names in the format required by decsg: a char matrix, ONE NAME PER COLUMN ---
%   names(:)  -> puts the strings in a column
%   char(...) -> converts them to a 16x2 char matrix (one name per row,
%                shorter names padded with a space)
%   '         -> transposes it to 2x16 (one name per column):
%                   M N P E P E P E ...
%                       1 1 2 2 3 3 ...
ns = char(names(:))';
 
% --- Set formula: how to combine the circles ---
%   +  union         -  difference (e.g. 'R1-C1' = rectangle with a hole)
%   *  intersection
% strjoin joins all the names with "+":  'M+N+P1+E1+P2+E2+...+P7+E7'
% We use the UNION because we want to KEEP every circle boundary: each one
% separates two different tissues. decsg keeps these internal boundaries and
% so divides the domain into many regions. (With a difference, the
% fascicles would become holes instead.)
sf = char(strjoin(names, "+"));
 
% --- Decomposition ---
% decsg (DEComposed Constructive Solid Geometry) splits the plane into the
% smallest regions delimited by the circles: these are the FACES of the geometry.
% Here there are 16 faces:
%   1 ring of surrounding medium + 1 region of epineurium (nerve minus fascicles)
%   + for each fascicle: 1 perineurium ring + 1 endoneurium disk  (7 x 2 = 14)
%
% Outputs:
%   dl = decomposed geometry: the description of the BOUNDARIES. Each column is
%        an arc of circle, with its coordinates and the faces on its left and
%        right. We do not need to read it: it is the input of fegeometry (section 5).
%   bt = Boolean table: bt(f, s) = 1 if face f lies INSIDE circle s.
%        One row per face, one column per circle (same order as gd). Examples:
%                                     M  N  P1 E1 P2 E2 ...
%          surrounding medium face    1  0  0  0  0  0
%          epineurium face            1  1  0  0  0  0
%          perineurium of fascicle 2  1  1  0  0  1  0
%          endoneurium of fascicle 2  1  1  0  0  1  1
%        The face NUMBERS are chosen by decsg, not by us: bt is what lets us
%        find out which tissue each face is (section 6).
[dl, bt] = decsg(gd, sf, ns);
 
 
%% ---------- 5. CREATE THE FINITE ELEMENT MODEL ----------
% femodel creates the model object that will contain everything:
% geometry, material properties, mesh (and, in Part 2, boundary conditions,
% sources and results). All of these are PROPERTIES of the variable "model".
%
% AnalysisType chooses the physics. "dcConduction" = current flowing through
% conductive materials (Ohm's law). Other options: "electrostatic",
% "thermalSteady", "structuralStatic", ...
%
% fegeometry(dl) converts the output of decsg into a geometry object.
%
% Here we only SET UP the model: no boundary conditions, no solution.
model  = femodel(AnalysisType="dcConduction", Geometry=fegeometry(dl));
nFaces = model.Geometry.NumFaces;
fprintf('\nGeometry: %d faces (regions), %d edges\n', nFaces, model.Geometry.NumEdges);
 
 
%% ---------- 6. TISSUE OF EACH GEOMETRY FACE ----------
% We need to know which tissue each face is, to give it the right material.
% We read it from bt (columns: 1 = M, 2 = N, then P1, E1, P2, E2, ...):
%   - inside an E circle               -> endoneurium
%   - inside a P circle but not an E   -> perineurium
%   - inside N but in no fascicle      -> epineurium
%   - otherwise (inside M only)        -> surrounding medium
% The ORDER of the checks matters: an endoneurium face is also inside its
% P circle and inside N, so we must test the innermost circle (E) first.
assert(size(bt,1) == nFaces, 'Unexpected number of faces.');   % stops with an error if false
colPeri = 3:2:nShapes;        % columns of the P circles: 3, 5, 7, ...
colEndo = 4:2:nShapes;        % columns of the E circles: 4, 6, 8, ...
 
faceTissue = zeros(nFaces, 1);     % tissue number of each face (1...4)
faceFasc   = zeros(nFaces, 1);     % fascicle number (0 = not in a fascicle)
for f = 1:nFaces
    if any(bt(f, colEndo))                       % any(): true if at least one element is 1
        faceTissue(f) = 1;   faceFasc(f) = find(bt(f, colEndo));   % find(): position of the 1
    elseif any(bt(f, colPeri))
        faceTissue(f) = 2;   faceFasc(f) = find(bt(f, colPeri));
    elseif bt(f, 2)
        faceTissue(f) = 3;
    else
        faceTissue(f) = 4;
    end
end
 
% Print the face IDs of each tissue (compare with the face-label plot, section 10)
for t = 1:4
    fprintf('   %-20s faces: %s\n', tissueName(t), mat2str(find(faceTissue == t)'));
end
 
 
%% ---------- 7. ASSIGN MATERIAL PROPERTIES ----------
% materialProperties creates an object with the physical properties of a
% material; model.MaterialProperties(faceIDs) = ... assigns it to those faces.
% This is the "material" information that every element of the mesh will
% carry (each element inherits the properties of the face it belongs to).
%
% Why only sigmaT? a good quesiton for the students
% 
% In a 2D cross-section the current can only flow in the
% x-y plane, i.e. ACROSS the nerve, so the transverse conductivity is the
% relevant one. sigmaL (along the nerve, z axis) would be needed in a 3D model.
% (Note: the toolbox accepts one isotropic value per region for
%  ElectricalConductivity.)
for t = 1:4
    faceIDs = find(faceTissue == t)';            % faces of tissue t
    model.MaterialProperties(faceIDs) = ...
        materialProperties(ElectricalConductivity = tissues.sigmaT(t));
end
fprintf('Material properties assigned to %d faces\n', nFaces);
 
 
%% ---------- 8. GENERATE THE MESH ----------
% generateMesh splits every face into triangles. Options used:
%   Hmax           = maximum element size (used where nothing else is specified,
%                    i.e. in the surrounding medium)
%   Hmin           = minimum element size
%   Hgrad          = growth rate between neighboring elements
%   Hface          = target size for groups of faces:
%                    {faceIDs, size, faceIDs, size, ...}  -> local refinement
%   GeometricOrder = "linear": triangles with 3 nodes (the vertices)
%                    "quadratic" (default): 6 nodes (vertices + edge midpoints),
%                    more accurate but harder to plot by hand
% ======================= EXERCISE 2 (part a) =======================
% Complete Hface: it must contain, for the endoneurium, the perineurium and
% the epineurium, the list of face IDs followed by the target element size
% (the sizes are defined in section 3).
% The first line (endoneurium) is given as an example.
% Hints:
%   - faceTissue (section 6) gives the tissue number of every face:
%       1 = endoneurium, 2 = perineurium, 3 = epineurium, 4 = surrounding medium
%   - find(faceTissue == t)' returns the IDs of the faces of tissue t (as a row)
%   - the surrounding medium is NOT listed: it will use Hmax
Hface = { find(faceTissue == 1)', hEndo, ...      % endoneurium (example)
          find(faceTissue == 2)', hPeri, ...      % perineurium
          find(faceTissue == 3)', hEpi };         % epineurium
 
% Automatic check: Hface must contain 3 pairs {faceIDs, size} = 6 entries
assert(numel(Hface) == 6, ...
       'EXERCISE 2a: Hface must contain the face IDs and the size for 3 tissues.');
 
% ======================= EXERCISE 2 (part b) =======================
% The line below generates a mesh with the DEFAULT options: it runs, but the
% elements are NOT the sizes we want (look at the table printed in section 9,
% and at the mesh plots: is the perineurium described well?).
% Add the options to generateMesh, using the parameters of section 3:
%   Hmax           -> maximum element size (used in the surrounding medium)
%   Hmin           -> minimum element size: use half of the perineurium size
%   Hgrad          -> growth rate between neighboring elements
%   Hface          -> the cell array you completed above
%   GeometricOrder -> "linear" (triangles with 3 nodes)
% Syntax example:  model = generateMesh(model, Hmax = 0.5, Hgrad = 1.5);
model = generateMesh(model, ...
                     Hmax           = hMedium, ...    % used where Hface says nothing (medium)
                     Hmin           = hPeri/2, ...    % half of the perineurium size
                     Hgrad          = hGrad, ...      % smooth growth fine -> coarse
                     Hface          = Hface, ...      % local sizes per tissue
                     GeometricOrder = "linear");      % 3-node triangles
% ====================================================================
 
% How the mesh is stored in memory:
%   Nodes    = 2 x Nnodes matrix: column k = [x; y] coordinates of node k
%   Elements = 3 x Nelements matrix (CONNECTIVITY): column e = IDs of the
%              3 nodes that form triangle e. The element does not store
%              coordinates: it points to its nodes.
msh   = model.Mesh;
nodes = msh.Nodes;
elems = msh.Elements;
nElem = size(elems, 2);       % size(A,2) = number of columns
fprintf('\nMesh: %d nodes, %d triangles\n', size(nodes,2), nElem);
fprintf('Example - element 1 is made of nodes %s\n', mat2str(elems(:,1)'));
 
 
%% ---------- 9. LINK EACH MESH ELEMENT TO ITS TISSUE (+ CHECK) ----------
% Every element belongs to exactly one geometry face, so it inherits that
% face's tissue. findElements(msh, "region", "Face", f) returns the IDs of
% all the elements inside face f.
%
% As a safety CHECK, for one element per face we recompute the tissue from
% its POSITION (distance of its centroid from the circle centers) and compare
% it with the tissue found from bt. If they differ, something is wrong.
elemTissue = zeros(1, nElem);      % tissue number of every element
nWrong = 0;
for f = 1:nFaces
    el = findElements(msh, "region", "Face", f);
    elemTissue(el) = faceTissue(f);
 
    c = mean(nodes(:, elems(:, el(1))), 2);                  % centroid = mean of its 3 nodes
    dFasc = hypot(c(1) - fasc(:,1), c(2) - fasc(:,2));       % distance from each fascicle center
                                                             % (hypot(a,b) = sqrt(a^2 + b^2))
    if     any(dFasc <= fasc(:,3)),          tCheck = 1;     % inside an endoneurium
    elseif any(dFasc <= fasc(:,3) + tPeri),  tCheck = 2;     % inside a perineurium ring
    elseif hypot(c(1), c(2)) <= Rnerve,      tCheck = 3;     % inside the nerve
    else,                                    tCheck = 4;     % outside the nerve
    end
    if tCheck ~= faceTissue(f), nWrong = nWrong + 1; end     % ~= means "not equal"
end
if nWrong == 0
    fprintf('Tissue check: OK\n');
else
    warning('%d faces out of %d have an unexpected tissue.', nWrong, nFaces);
end
 
% --- Element size and quality ---
% Area of each triangle from the coordinates of its 3 nodes (vectorized:
% all triangles at once, no loop). i1, i2, i3 = first, second, third node
% of every element.
x = nodes(1,:);   y = nodes(2,:);
i1 = elems(1,:);  i2 = elems(2,:);  i3 = elems(3,:);
elemArea = 0.5 * abs((x(i2)-x(i1)).*(y(i3)-y(i1)) - (x(i3)-x(i1)).*(y(i2)-y(i1)));
elemSize = sqrt(4*elemArea/sqrt(3));   % edge of the equilateral triangle with the same area
quality  = meshQuality(msh);           % shape quality: 1 = perfect equilateral, 0 = degenerate
                                       % (very flat elements give inaccurate results)
 
% Summary: is the mesh as fine as we asked, in every tissue?
hTarget = [hEndo; hPeri; hEpi; hMedium];
fprintf('\n%-20s %10s %14s %14s\n', 'Tissue', 'Elements', 'Target h [mm]', 'Mean h [mm]');
for t = 1:4
    m = elemTissue == t;                         % logical mask: true for elements of tissue t
    fprintf('%-20s %10d %14.3f %14.3f\n', tissueName(t), nnz(m), hTarget(t), mean(elemSize(m)));
end
fprintf('Mesh quality: mean %.2f, min %.2f\n', mean(quality), min(quality));
 
 
%% ---------- 10. PLOT: GEOMETRY ----------
% Three views: (a) the anatomy drawn with colored circles, (b) the geometry
% as the toolbox sees it, with the face numbers, (c) a 3D sketch of the nerve.
th = linspace(0, 2*pi, 200);       % 200 angles from 0 to 2*pi, used to draw circles:
                                   % a circle is x = xc + r*cos(th), y = yc + r*sin(th)
 
figure(Color="w", Position=[40 80 1500 480], Name="Geometry");
tiledlayout(1, 3, TileSpacing="compact");   % 1 row x 3 columns of plots
 
% 10a. Anatomy of the cross-section
% fill(x, y, color) draws a filled polygon. Circles are drawn from the
% outside to the inside, so that the inner ones are drawn on top.
nexttile
fill(Rnerve*1.25*cos(th), Rnerve*1.25*sin(th), tissueColor(4,:), EdgeColor="none"); hold on
fill(Rnerve*cos(th), Rnerve*sin(th), tissueColor(3,:), EdgeColor=[0.6 0.45 0.2], LineWidth=1.5);
for i = 1:nFasc
    rOut = fasc(i,3) + tPeri;
    fill(fasc(i,1) + rOut*cos(th),      fasc(i,2) + rOut*sin(th),      tissueColor(2,:), EdgeColor="none");
    fill(fasc(i,1) + fasc(i,3)*cos(th), fasc(i,2) + fasc(i,3)*sin(th), tissueColor(1,:), EdgeColor="none");
    text(fasc(i,1), fasc(i,2), num2str(i), HorizontalAlignment="center", FontWeight="bold");
end
axis equal; box on                 % axis equal: same scale on x and y (circles look round)
xlim(1.2*Rnerve*[-1 1]);  ylim(1.2*Rnerve*[-1 1]);
xlabel('x [mm]'); ylabel('y [mm]');
title('Nerve cross-section (fascicle numbers)');
 
% 10b. Geometry as seen by the PDE Toolbox, with face numbers.
% These are the IDs used to assign materials and mesh sizes.
nexttile
pdegplot(model.Geometry, FaceLabels="on");
% pdegplot draws its lines with clipping OFF: when we zoom on the nerve, the
% outer medium circle would be drawn outside the axes, over the other plots.
% Turning clipping ON cuts everything at the axes borders.
set(findall(gca, '-property', 'Clipping'), 'Clipping', 'on');
axis equal
xlim(1.2*Rnerve*[-1 1]);  ylim(1.2*Rnerve*[-1 1]);
xlabel('x [mm]'); ylabel('y [mm]');
title('Geometry faces (IDs used by the model)');
 
% 10c. 3D view of a nerve segment (drawing only, not part of the model).
% cylinder(1, 100) returns the surface of a cylinder with radius 1 and
% height 1 (100 points around it); we scale and shift it for every tissue.
nexttile
[Xc, Yc, Zc] = cylinder(1, 100);
Zc = Lnerve * Zc;                  % stretch the height to Lnerve
surf(Rnerve*Xc, Rnerve*Yc, Zc, FaceColor=tissueColor(3,:), FaceAlpha=0.25, EdgeColor="none"); hold on
for i = 1:nFasc
    rOut = fasc(i,3) + tPeri;
    surf(fasc(i,1) + rOut*Xc, fasc(i,2) + rOut*Yc, Zc, ...
         FaceColor=tissueColor(2,:), FaceAlpha=0.7, EdgeColor="none");
end
% End cap = the cross-section at z = Lnerve, drawn with fill3 (fill in 3D).
% Small offsets in z (0.01, 0.02) keep the layers from overlapping visually.
zTop = Lnerve * ones(size(th));
fill3(Rnerve*cos(th), Rnerve*sin(th), zTop, tissueColor(3,:), EdgeColor=[0.6 0.45 0.2]);
for i = 1:nFasc
    rOut = fasc(i,3) + tPeri;
    fill3(fasc(i,1) + rOut*cos(th),      fasc(i,2) + rOut*sin(th),      zTop + 0.01, tissueColor(2,:), EdgeColor="none");
    fill3(fasc(i,1) + fasc(i,3)*cos(th), fasc(i,2) + fasc(i,3)*sin(th), zTop + 0.02, tissueColor(1,:), EdgeColor="none");
end
axis equal; grid on; view(-35, 20);   % view(azimuth, elevation) sets the 3D viewpoint
xlabel('x [mm]'); ylabel('y [mm]'); zlabel('z [mm]');
title(sprintf('3D view of a %g mm nerve segment', Lnerve));
 
 
%% ---------- 11. PLOT: MESH COLORED BY TISSUE ----------
% patch draws all the triangles at once:
%   Faces           = which nodes form each triangle (the connectivity, transposed)
%   Vertices        = node coordinates (transposed)
%   FaceVertexCData = one color per triangle: the color of its tissue
%   FaceColor="flat" = each triangle has a single uniform color
% Three views of the same mesh: whole domain, whole nerve, zoom on fascicle 2.
% Note how the elements are small in the perineurium and large far away.
viewX  = {Rmedium*[-1 1], 1.1*Rnerve*[-1 1], fasc(2,1) + 0.5*[-1 1]};   % x limits of each view
viewY  = {Rmedium*[-1 1], 1.1*Rnerve*[-1 1], fasc(2,2) + 0.5*[-1 1]};   % y limits of each view
viewTitle = ["Whole domain", "Nerve", "Zoom on fascicle 2"];
 
figure(Color="w", Position=[60 60 1500 500], Name="Mesh");
tiledlayout(1, 3, TileSpacing="compact");
for v = 1:3
    nexttile
    patch(Faces=elems', Vertices=nodes', FaceVertexCData=tissueColor(elemTissue,:), ...
          FaceColor="flat", EdgeColor=[0.25 0.25 0.25], LineWidth=0.2);
    axis equal; box on
    xlim(viewX{v});  ylim(viewY{v});
    xlabel('x [mm]'); ylabel('y [mm]');
    title(viewTitle(v));
end
% Legend: one invisible patch (NaN coordinates) per tissue, used only for
% its color in the legend.
hold on
hLeg = gobjects(4,1);
for t = 1:4
    hLeg(t) = patch(NaN, NaN, tissueColor(t,:));
end
legend(hLeg, tissueName, Location="southoutside", Orientation="horizontal");
 
 
%% ---------- 12. PLOT: MESH STATISTICS ----------
% Checking the mesh before solving is good practice:
%   (a) how many elements in each tissue
%   (b) is the element size close to what we asked? (box = spread of sizes,
%       red x = target size from section 3)
%   (c) are the elements well shaped? (quality close to 1 is good)
figure(Color="w", Position=[80 80 1400 420], Name="Mesh statistics");
tiledlayout(1, 3, TileSpacing="compact");
tissueCat = categorical(tissueName(elemTissue(:)), tissueName);   % tissue name of every element,
                                                                  % keeping the tissue order
nexttile
% histcounts with edges 0.5, 1.5, ..., 4.5 counts how many elements have tissue 1, 2, 3, 4
bar(categorical(tissueName, tissueName), histcounts(elemTissue, 0.5:1:4.5));
ylabel('Number of elements'); title('Elements per tissue');
 
nexttile
boxchart(tissueCat, elemSize(:)); hold on
plot(categorical(tissueName, tissueName), hTarget, 'rx', MarkerSize=12, LineWidth=2);
set(gca, YScale="log");            % log scale: sizes span more than one order of magnitude
ylabel('Element size [mm]'); title('Element size (red x = target)');
 
nexttile
histogram(quality, 30);            % 30 bins
xlabel('Quality (1 = equilateral)'); ylabel('Number of elements');
title('Mesh quality');
 
 
%% ---------- 13. PLOT: ELECTRICAL PROPERTY MAPS ON THE MESH ----------
% Each element is colored by the property of its tissue: this shows the
% "material information" carried by every element of the mesh.
% A log10 scale is used because the values span several orders of magnitude
% (the perineurium is ~100 times less conductive than the other tissues).
% Note: the longitudinal map shows what each element would carry in a 3D
% model; in this 2D model only sigmaT is used.
propValue = [tissues.sigmaT, tissues.sigmaL, tissues.rhoT];   % 4 x 3: tissues x properties
propTitle = ["Transverse conductivity \sigma_T [S/m]", ...
             "Longitudinal conductivity \sigma_L [S/m]", ...
             "Transverse resistivity \rho_T [\Omega m]"];
 
figure(Color="w", Position=[100 60 1500 480], Name="Electrical properties");
tiledlayout(1, 3, TileSpacing="compact");
for k = 1:3
    nexttile
    logVal = log10(propValue(elemTissue, k));          % one value per element
    patch(Faces=elems', Vertices=nodes', FaceVertexCData=logVal, ...
          FaceColor="flat", EdgeColor="none");
    axis equal; box on
    xlim(1.1*Rnerve*[-1 1]);  ylim(1.1*Rnerve*[-1 1]);
    xlabel('x [mm]'); ylabel('y [mm]');
    colormap(gca, parula);
    cb = colorbar;
    % Colorbar ticks at whole powers of 10, labeled as 10^n
    tk = floor(min(logVal)):ceil(max(logVal));
    clim([tk(1) tk(end)]);
    cb.Ticks = tk;
    cb.TickLabels = "10^{" + tk + "}";
    title(propTitle(k));
end
 
 
%% ---------- 14. PLOT: TISSUE RESISTIVITY ----------
% Resistivity rho = 1/sigma says how much a tissue OPPOSES the current flow
% (high rho = poor conductor). Comparing the transverse and longitudinal
% values shows the ANISOTROPY: in the endoneurium and in the surrounding
% medium the current flows more easily along the nerve (lower rhoL).
% A log scale is needed because the perineurium is ~100 times more resistive
% than the other tissues.
figure(Color="w", Position=[120 60 700 450], Name="Tissue resistivity");
 
% ======================= EXERCISE 3 =======================
% Make a BAR CHART that compares, for each tissue, the transverse (rhoT) and
% the longitudinal (rhoL) resistivity: two bars side by side for each tissue.
% The function to create this plot is "bar".
% Steps:
%   1. x values: the tissue names. Use  categorical(tissueName, tissueName)
%      (categorical turns text into categories; the second input keeps
%       the tissues in our order instead of alphabetical order)
%   2. y values: a 4 x 2 matrix, one column for rhoT and one for rhoL
%      (with a matrix, bar draws one group of bars for each row)
%   3. logarithmic y axis:  set(gca, YScale="log")
%      Question: why is a log scale needed here?
%   4. add: a grid, the y-axis label with units (resistivity in Ohm*m;
%      tip: '\Omega' writes the Greek letter), a legend
%      ('Transverse', 'Longitudinal') and a title
% Useful functions: bar, categorical, set, grid, ylabel, legend, title
% (type "doc bar" in the Command Window to see examples)
 
% Answer to the question: the resistivities span more than 3 orders of
% magnitude (perineurium ~1150 Ohm*m, the other tissues ~2-12 Ohm*m). On a
% linear axis only the perineurium bars would be visible.
tissueCatX = categorical(tissueName, tissueName);       % x: tissues, in our order
bar(tissueCatX, [tissues.rhoT, tissues.rhoL]);          % 4x2 -> 2 bars per tissue
set(gca, YScale="log");
grid on
ylabel('Resistivity \rho [\Omega m]');
legend('Transverse', 'Longitudinal', Location="northeast");
title('Tissue resistivity: transverse vs longitudinal');

 
% ==========================================================
 
 
%% ---------- 15. FIGURE APPEARANCE: DARK TEXT ON WHITE BACKGROUND ----------
% When MATLAB uses a DARK theme, figures get light-gray text and dark axes,
% which are hard to read on a projector or in a report. Here we:
%   1. switch every figure to the light theme (MATLAB R2025a and later)
%   2. force dark colors on axes, titles, labels, legends and colorbars
%      (works in any MATLAB version)
% findall(fig, 'Type', ...) finds all the objects of a given type in a figure,
% including titles and axis labels (which are "text" objects).
allFigs = findall(0, 'Type', 'figure');      % 0 = the MATLAB root: all open figures
for fig = allFigs'
    if exist('theme', 'file')
        theme(fig, "light");
    end
    set(fig, 'Color', 'w');
    set(findall(fig, 'Type', 'axes'),     'Color', 'w', ...
        'XColor', 'k', 'YColor', 'k', 'ZColor', 'k', 'GridColor', [0.3 0.3 0.3]);
    set(findall(fig, 'Type', 'text'),     'Color', 'k');   % titles, axis labels, numbers
    set(findall(fig, 'Type', 'legend'),   'TextColor', 'k', 'Color', 'w', 'EdgeColor', 'k');
    set(findall(fig, 'Type', 'colorbar'), 'Color', 'k');
end
 
fprintf('\nDone: geometry, material properties and mesh are ready.\n');
fprintf('Next steps (Part 2): boundary conditions, sources and solve.\n');


%Geometry of the node with electrical site in it. Geometry important and a
%simulation and plot of it.
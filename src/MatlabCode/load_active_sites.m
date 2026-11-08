%% ACTIVE SITES OF THE NERVE INTERFACE
% This script loads the matrix "sites", which contains the position, the
% orientation and the size of the 48 active sites of the nerve interface.
% Use these data to insert the active sites in your 3D nerve model.
%
% ---------------------------------------------------------------------------
% REFERENCE FRAME AND UNITS
% ---------------------------------------------------------------------------
%   z axis   = nerve axis; origin at the centre of the device
%   x-y      = cross-section plane of the nerve
%   r, theta = cylindrical coordinates in the x-y plane:
%              x = r*cosd(theta),  y = r*sind(theta)
%              theta measured from the x axis, counterclockwise
%   Units: mm for lengths, degrees for angles
%
% ---------------------------------------------------------------------------
% THE DEVICE
% ---------------------------------------------------------------------------
%   The device has 8 SECTIONS along the nerve, numbered -4, -3, -2, -1
%   (z < 0) and +1, +2, +3, +4 (z > 0). The sections are 2 mm apart and
%   sections -1 and +1 are at 4.3 mm from z = 0. Section +k is the mirror
%   image of section -k with respect to the plane z = 0.
%
%   Each section contains 6 active sites:
%     - 4 INTRANEURAL sites, on a shank that crosses the whole nerve along a
%       diameter. Their distance from the nerve axis is fixed by the device
%       (0.3 and 0.9 mm): it does NOT depend on the nerve radius.
%     - 2 CUFF sites, on the surface of the nerve. Their distance from the
%       axis DEPENDS on the radius of your nerve.
%
% ---------------------------------------------------------------------------
% STRUCTURE OF THE MATRIX "sites"   (48 rows x 10 columns)
% ---------------------------------------------------------------------------
%   One ROW for each active site:
%
%    1  Electrode   electrode number (1...48)
%    2  Type        1 = intraneural site,  2 = cuff site
%    3  Section     section number: -4...-1 (z<0) or +1...+4 (z>0)
%    4  z           position of the centre along the nerve axis [mm]
%    5  r           distance of the centre from the nerve axis [mm]
%                   NaN for the cuff sites: they lie ON the nerve surface,
%                   so for them use  r = R + thickness/2  (R = your nerve radius),
%                   so that the face of the site rests on the nerve surface
%    6  theta       angular position of the centre [deg]
%    7  phi         orientation of the site: direction of the NORMAL to its
%                   square face, i.e. the direction of the thickness [deg].
%                   The normal always lies in the x-y plane:
%                       n = [cosd(phi), sind(phi), 0]
%                   - intraneural sites: the normal is perpendicular to the
%                     shank (the face lies along the diameter and the nerve
%                     axis: the thin side slides between the fibres)
%                   - cuff sites: the normal points towards the nerve axis
%    8  thickness   thickness of the site, along the normal [mm]
%    9  sideZ       side of the square face parallel to the nerve axis [mm]
%   10  sideXY      side of the square face lying in the x-y plane [mm],
%                   along the direction  w = [-sind(phi), cosd(phi), 0]
%
%   Each site is therefore a box with:
%     centre      c = [r*cosd(theta), r*sind(theta), z]
%     thickness   along n = [cosd(phi), sind(phi), 0]
%     sideZ       along the z axis
%     sideXY      along w = [-sind(phi), cosd(phi), 0]
%
%   Example: all the intraneural sites of section -2
%       sites(sites(:,2) == 1 & sites(:,3) == -2, :)
% ---------------------------------------------------------------------------

load("active_sites.mat", "sites");


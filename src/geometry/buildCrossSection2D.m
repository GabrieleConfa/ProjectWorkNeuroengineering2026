function xsec = buildCrossSection2D(patient, params)
%BUILDCROSSSECTION2D 2D geometry of the nerve cross-section with the tissue of each face.
%   xsec = buildCrossSection2D(patient, params)
%
%   Owner: WP1. Status: SKELETON — start from sections 4-6 of the tutor's
%   src/MatlabCode/PW_lesson_5_11_students.m (decsg + face -> tissue map), in metres.
%
%   Inputs
%     patient  struct from makePatient
%     params   struct from default_params() (uses params.domain.outerRadiusFactor)
%
%   Output: struct
%     dl            decomposed geometry from decsg [m] (input of fegeometry / extrude)
%     faceTissue    (nFaces,1) 1 endo, 2 peri, 3 epi, 4 medium
%     faceFascicle  (nFaces,1) fascicle index, 0 outside the fascicles
%     outerRadius_m radius of the surrounding medium [m]
%
%   Test to write: the tutor's position check (section 9) as a unit test.

arguments
    patient (1,1) struct
    params (1,1) struct
end
error('project:notImplemented', 'buildCrossSection2D: not implemented yet (WP1).');
end

function xsec = buildCrossSection2D(patient, params)
%BUILDCROSSSECTION2D 2D geometry of the nerve cross-section with the tissue of each face.
%
%   NOT IMPLEMENTED (final status, 2026-11-28). Skeleton from the first work plan (docs/work_plan.md),
%   never filled in: the project did this job elsewhere. Calling it stops with an error. Instead:
%   the 2D model builds its section inside src/main/main_pipeline.m; the 3D model builds the patient with
%   makePatient / gmshAnatomy and the geometry in src/main/build3DModel.m.
%   Kept only as a record of the planned interface; not used by any pipeline.
%
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

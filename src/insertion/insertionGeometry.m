function q = insertionGeometry(patient, p)
%INSERTIONGEOMETRY Anatomy of the main 3D model in the input format of the insertion code.
%   q = insertionGeometry(patient, p)
%
%   Purpose
%     The lead-insertion Monte Carlo (src/insertion/simulateInsertion.m, randomFascicleShapes.m) comes
%     from Marco's branch cuff-oval-viewer-convergence, where the anatomy is described by his
%     quicklook parameter struct q. This function fills the few fields that code reads from OUR
%     virtual patient (makePatient, config/pw6_3d_params.m), so the same insertion code runs on the
%     patient of main_pipeline_3d (merged 2026-11-26).
%
%   Inputs
%     patient  struct from makePatient (oval nerve: nerveSemiAxes_m, nerveRotation_deg; circular
%              fascicles: fascicles.x_m, y_m, rEndo_m = endoneurium radius)
%     p        parameters (config/pw6_3d_params.m)
%
%   Output: struct q
%     q.geometry.fascicles_m        [x y a b angle] per fascicle (circles: a = b = rEndo, angle 0) [m]
%     q.geometry.nerveSemiAxes_m    [a b] of the oval nerve [m]
%     q.geometry.nerveRotation_deg  rotation of the major axis [deg]
%     q.geometry.radius_m           equivalent nerve radius [m] (loadActiveSites)
%     q.geometry.tPeri_m            perineurium thickness [m]: the tip touches a fascicle when its
%                                   centre is within tip radius + tPeri of the endoneurium outline
%     q.electrode.sitesFile         the provided contact layout
%     q.electrode.reverseEntrySections  stations whose lead enters from the other side: none here
%                                   [assumed]; the entry side is not known (docs/project_brief.md)

arguments
    patient (1,1) struct
    p (1,1) struct
end
F = patient.fascicles;
q = struct();
q.geometry.fascicles_m = [F.x_m, F.y_m, F.rEndo_m, F.rEndo_m, zeros(height(F), 1)];
q.geometry.nerveSemiAxes_m = patient.nerveSemiAxes_m;
q.geometry.nerveRotation_deg = patient.nerveRotation_deg;
q.geometry.radius_m = patient.nerveRadius_m;
q.geometry.tPeri_m = p.perineurium.realThickness_m;
q.electrode.sitesFile = p.electrode.sitesFile;
q.electrode.reverseEntrySections = [];
end

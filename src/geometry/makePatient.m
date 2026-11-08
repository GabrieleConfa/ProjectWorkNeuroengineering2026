function patient = makePatient(id, seed, params)
%MAKEPATIENT Create a virtual patient: nerve radius and fascicles of the cross-section.
%   patient = makePatient(id, seed, params)
%
%   Owner: WP1. Status: SKELETON — contract agreed in docs/architecture.md §4.
%
%   Inputs
%     id      string, e.g. "P00"
%     seed    integer for rng(seed); NaN -> the tutor's fixed section
%             (params.geometry.tutorSection), i.e. patient P00
%     params  struct from default_params()
%
%   Output: struct
%     id, seed, level
%     nerveRadius_m   nerve radius [m]
%     fascicles       table: x_m, y_m, rEndo_m, tPeri_m  [m]
%
%   Sprint 1: implement the seed = NaN case (P00). Sprint 3: random generator
%   from the Delgado-Martinez 2016 statistics (no overlaps, inside the nerve).
%   Save with: save(fullfile('data','patients', id + ".mat"), 'patient')

arguments
    id (1,1) string
    seed (1,1) double
    params (1,1) struct
end
error('project:notImplemented', 'makePatient: not implemented yet (WP1).');
end

function run_personal_3d() % Run the full personal case from the MATLAB command window.
%% 1. Select the personal anatomy for this run
previous = getenv('NERVE_CASE'); % Remember the previous patient selection.
restore = onCleanup(@()setenv('NERVE_CASE',previous)); % Restore the selection on completion or error.
setenv('NERVE_CASE','personal'); % Select the 60 mm geometry and the routed pads.
%% 2. Generate tetrahedra, verify all cuff interfaces and solve the 48 fields
runStage('main_pipeline_3d'); % Mesh the volume and solve all contacts.
root = fileparts(fileparts(mfilename('fullpath'))); % Locate the saved result independently of the caller directory.
%% 3. Validate current balance and render the personal mesh and station sheets
check_personal_field(fullfile(root,'results','main3d','Cuff60_3d.mat')); % Audit the actual solved fields and locations.
plot_solved_stations(fullfile(root,'results','main3d','Cuff60_3d.mat')); % Generate all eight field sheets with the actual lead overlays.
%% 4. Export the same solved device for the interactive viewer
export_personal_preview(); % Refresh the geometry used by personal viewer links.
runStage('export_viewer_data'); % Export geometry and fields from this exact solution.
end % Finish the personal simulation.
function runStage(name) % Isolate scripts that clear their workspace.
run(which(name)); % Execute the stage with project paths already set by setup_paths.
end % Keep the parent cleanup and selection intact.

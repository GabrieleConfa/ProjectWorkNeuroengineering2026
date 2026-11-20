function [fig, model] = run_complete_preview(mode, anatomyParams) % Generate the existing 2D figures and the new 3D preview.
arguments % Accept the same editable inputs as the 3D entry point.
    mode (1,1) string {mustBeMember(mode,["deflect","penetrate"])} = "deflect" % Choose the insertion scenario.
    anatomyParams (1,1) struct = quicklook3DParams() % Read all provisional geometry and device inputs.
end % Finish the input declaration.
%% 1. Preserve the existing 2D graphics pipeline
runExisting2D(); % Run all original geometry, mesh, potential, field and fascicle plots.
fprintf("Historical 2D plots use the tutor P00 seven-fascicle anatomy; see results/2d_sites.\n"); % Label the existing outputs accurately.
%% 2. Generate the 3D anatomy and contact graphics
[fig,model] = run_3d_quicklook(mode,anatomyParams); % Build and draw the 3D configuration.
outDir = fullfile(fileparts(fileparts(mfilename('fullpath'))),'results','3d_preview'); % Locate the 3D output folder.
if ~exist(outDir,'dir'), mkdir(outDir); end % Create the output folder on the first run.
exportgraphics(fig,fullfile(outDir,"nerve_"+mode+".png"),Resolution=180); % Save the new 3D and end-on views.
savefig(fig,fullfile(outDir,"nerve_"+mode+".fig")); % Save an interactive MATLAB figure.
exportgraphics(model.graphics.atlas,fullfile(outDir,"all_stations_"+mode+".png"),Resolution=180); % Save all eight current anatomical sections.
exportgraphics(model.graphics.deviation,fullfile(outDir,"deflection_"+mode+".png"),Resolution=180); % Save the routed versus straight comparison.
exportgraphics(model.graphics.cuff,fullfile(outDir,"local_cuff_"+mode+".png"),Resolution=180); % Save the magnified cuff footprint and axial return.
end % Return the 3D figure and model after all graphics are produced.
function runExisting2D() % Isolate the original script's clearvars from the new entry point.
root = fileparts(fileparts(mfilename('fullpath'))); % Locate the repository root.
run(fullfile(root,'scripts','run_2d_sites_quicklook.m')); % Produce every figure from the existing 2D script.
end % Keep the new caller's parameters intact.

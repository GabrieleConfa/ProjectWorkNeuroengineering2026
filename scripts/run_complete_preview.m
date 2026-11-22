function [fig, model] = run_complete_preview(mode, anatomyParams) % Generate the existing 2D figures and the new 3D preview.
arguments % Accept the same editable inputs as the 3D entry point.
    mode (1,1) string {mustBeMember(mode,["deflect","penetrate"])} = "deflect" % Choose the insertion scenario.
    anatomyParams (1,1) struct = quicklook3DParams() % Read all provisional geometry and device inputs.
end % Finish the input declaration.
%% 1. Draw the current device layout
plot_personal_layout(anatomyParams); % Use the same anatomy as the 3D preview.
%% 2. Generate the 3D anatomy and contact graphics
[fig,model] = run_3d_quicklook(mode,anatomyParams); % Build and draw the 3D configuration.
outDir = fullfile(fileparts(fileparts(mfilename('fullpath'))),'results','3d_preview'); % Locate the 3D output folder.
if ~exist(outDir,'dir'), mkdir(outDir); end % Create the output folder on the first run.
exportgraphics(fig,fullfile(outDir,"nerve_"+mode+".png"),Resolution=180); % Save the new 3D and end-on views.
savefig(fig,fullfile(outDir,"nerve_"+mode+".fig")); % Save an interactive MATLAB figure.
exportgraphics(model.graphics.atlas,fullfile(outDir,"all_stations_"+mode+".png"),Resolution=180); % Save all eight current anatomical sections.
if mode=="deflect" % Keep the personal review link in sync on future geometry changes.
    copyfile(fullfile(outDir,"all_stations_deflect.png"),fullfile(outDir,"figura2_equilibrio.png")); % Reuse the exact exported atlas.
end % Preserve the distinct penetration scenario.
exportgraphics(model.graphics.deviation,fullfile(outDir,"deflection_"+mode+".png"),Resolution=180); % Save the routed versus straight comparison.
exportgraphics(model.graphics.cuff,fullfile(outDir,"local_cuff_"+mode+".png"),Resolution=180); % Save the magnified cuff footprint and axial return.
model.graphics.functions = plotFascicleMap(anatomyParams); % Draw dimensions and available functional assignments.
exportgraphics(model.graphics.functions,fullfile(outDir,"fascicle_map.png"),Resolution=180); % Save the labelled anatomy beside the eight sections.
end % Return the 3D figure and model after all graphics are produced.

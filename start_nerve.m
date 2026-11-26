function start_nerve(action) % Press Run once to calculate and display the personal nerve model.
%% 1. Set up this checkout
if nargin==0,action="run";end % Run the complete model when launched from the MATLAB editor.
root=fileparts(mfilename('fullpath')); addpath(root); setup_paths(); % Resolve every source relative to this file.
previousVisibility=get(groot,'DefaultFigureVisible'); set(groot,'DefaultFigureVisible','off'); restoreVisibility=onCleanup(@()set(groot,'DefaultFigureVisible',previousVisibility)); %#ok<NASGU> % Draw figures off screen: they are saved as PNG files and shown in the browser gallery (this also keeps unattended runs reliable).
%% 2. Run the personal workflow or reopen existing results
switch string(action) % Shortcuts open a saved result without repeating the numerical solve.
    case "run"
        q=quicklook3DParams(); % Load the common anatomy, material and closure settings.
        fprintf('1/5  Anatomy, routed leads and figures\n'); run_complete_preview("deflect",q); % Draw the same section used in the FEM.
        fprintf('2/5  Connected sections, tetrahedral mesh and 48 electrical fields\n'); run_personal_3d(); % Solve the entire personal case.
        fprintf('3/5  Mesh convergence study (4 meshes, 3 pads)\n'); run_convergence_study(); % Check how much the results still change when the mesh is refined.
        fprintf('4/5  Random fascicle shapes and lead insertion (Monte Carlo)\n'); run_insertion_study(); % Scenario study with irregular fascicles.
        fprintf('5/5  Opening the gallery, mesh and stimulation viewer in Chrome\n'); % Report the final display step.
        pages=["figure","mesh","personal"]; % Open the field viewer last.
    case "convergence"
        run_convergence_study(); pages="figure"; % Repeat only the mesh study and show it in the gallery.
    case "insertion"
        run_insertion_study(); pages="figure"; % Repeat only the insertion study and show it in the gallery.
    case "geometry"
        run_complete_preview(); pages="figure"; % Rebuild the anatomical views alone.
    case "open"
        pages=["figure","mesh","personal"]; % Reopen the three saved pages without repeating the simulation.
    case "viewer"
        pages="personal"; % Open the solved stimulation case.
    case "mesh"
        pages="mesh"; % Open the tetrahedral mesh.
    case "figures"
        pages="figure"; % Open the current figure gallery.
    otherwise
        error('Use start_nerve, or specify open, convergence, insertion, geometry, figures, mesh or viewer.'); % Report unsupported operations.
end % Finish the requested numerical work.
%% 3. Open the interactive results in Chrome using an absolute, URL-encoded local path
for page=pages % Show each requested page.
    path=fullfile(root,'viewer',page+".html"); % Locate the offline HTML file.
    if ~isfile(path),warning('Missing %s: run start_nerve first.',path);continue;end % Report a missing page instead of opening a blank tab.
    localPath=replace(string(path),"\","/"); % Normalize Windows separators before URI encoding.
    parts=split(localPath,"/"); % Encode path components without escaping directory separators.
    for k=1:numel(parts),parts(k)=string(urlencode(char(parts(k))));end % Escape spaces, Unicode and reserved characters.
    encoded=replace(join(parts,"/"),"+","%20"); % Use path encoding for spaces.
    if ispc,encoded="/"+replace(encoded,"%3A",":");end % Retain the Windows drive prefix.
    url=char("file://"+encoded); % Build a local URL without requiring a Java runtime.
    fprintf('%s\n',url); % Keep the exact link available in the Command Window.
    openInChrome(url); % The viewer needs WebGL2 and local scripts, which Chrome handles reliably.
end % Finish displaying the results.
end
function openInChrome(url) % Open a local page in Google Chrome and fall back to the default browser.
if ismac,cmd=sprintf('open -a "Google Chrome" "%s"',url); % macOS application launcher.
elseif ispc,cmd=sprintf('start "" chrome "%s"',url); % Windows launcher.
else,cmd=sprintf('google-chrome "%s" >/dev/null 2>&1 &',url); % Linux launcher.
end % Finish choosing the platform command.
status=system(cmd); % Run the launcher (works also when MATLAB runs without the desktop).
if status~=0,warning('Chrome not found. Open this link in Chrome: %s',url);web(url,'-browser'); end % Fall back to the default browser and show the link.
end

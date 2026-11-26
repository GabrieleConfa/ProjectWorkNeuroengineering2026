function exportPersonalMeshView(o,folder) % Export the actual exterior faces of selected tetrahedra for the interactive mesh view.
%% 1. Locate each solid element in the saved MATLAB solution
X=o.mesh.nodes_m'*1000; % Convert positions to millimetres for the viewer.
T=double(o.mesh.tets'); % Store one solved tetrahedron per row.
material=o.cells.tissue(double(o.mesh.elemCell(:))); % Recover tissue and metal identities.
z=mean(reshape(X(T,3),size(T)),2); % Locate each element along the nerve.
y=mean(reshape(X(T,2),size(T)),2); % Choose whole elements for the longitudinal cutaway.
station=o.sites(o.sites.type=="cuff",:); % Read the physical electrode stations.
[sections,idx]=unique(station.section); % Select each axial location once.
centres=station.z_m(idx)*1000; % Use the same coordinates as the electrical solution.
data=struct('count',size(T,1),'sections',sections','centres',centres','views',{{}}); % Keep view labels and the full element count together.
%% 2. Collect visible faces of whole elements at each station and along the full nerve
[~,interfaces]=checkCuffMeshContact(X'/1000,T',material,o.cells.site(double(o.mesh.elemCell(:))),o.sites,o.params.personalGeometry); % Identify the actual shared face of every exterior pad.
contactFaces=sort(vertcat(interfaces{:}),2); % The epineurium owns the visible face when viewed from inside.
siteIds=o.cells.site(double(o.mesh.elemCell(:))); % Distinguish cuff metal from penetrating contacts.
external=ismember(siteIds,station.id); % Remove the corresponding cuff metal when the overlying tissue is cut away.
for k=0:numel(sections) % Include both the overall model and local station views.
    selected=true(size(z)); % Start with the complete solved length.
    if k>0, selected=abs(z-centres(k))<=0.85; end % Display a 1.7 mm window about this station.
    surfaces=struct(); % Keep each material independently selectable.
    for kind=[1 3 5 6] % Display endoneurium, epineurium, active metal and the PDMS lead bodies.
        inside=selected & material==kind; % Select tetrahedra of this material.
        name="t"+kind; % Use the solver tissue identifiers.
        omit=zeros(0,3); % Keep all faces of ordinary tissue boundaries.
        if kind==5, omit=contactFaces; end % Avoid a coincident metal face competing with the tissue.
        surfaces.(name)=surfaceData(X,T(inside,:),omit); % Export only faces exposed by this element selection.
        if kind==5,surfaces.cutMetal=surfaceData(X,T(inside & (~external | y<=0),:),omit);end % Remove external pads on the opened side.
        if kind==3,surfaces.cut=surfaceData(X,T(inside & y<=0,:));end % Reveal internal structures by removing half of the epineurial elements.
    end % Complete the anatomical layers.
    data.views{end+1}=surfaces; % Retain the full model first, then the stations in axial order.
end % Finish the mesh views.
%% 3. Save an offline JavaScript data file alongside the solver results
path=fullfile(folder,'Cuff60_mesh.js'); % Match the local gallery link.
fid=fopen(path+".tmp",'w'); assert(fid>=0,'Cannot write mesh view.'); % Stage output until the serialization completes.
guard=onCleanup(@()fclose(fid)); % Close the staged file on failure.
fprintf(fid,'window.nerveMesh=%s;',jsonencode(data)); % Write only surface buffers derived from actual tetrahedra.
clear guard; movefile(path+".tmp",path,'f'); % Publish the complete mesh in one operation.
%% 4. Render the first station for the gallery with the same cutaway and tissue colours
fig=figure(Visible="off",Color="white",Theme="light",Position=[50 50 880 570]); % Create a reproducible static 3D view.
ax=axes(fig); hold(ax,'on'); % Group the mesh layers on one set of physical axes.
keys=["cut","t1","cutMetal","t6"]; colours=[0.92 0.78 0.58;0.88 0.50 0.42;0.55 0.22 0.74;0.48 0.65 0.92]; % Match the tissue, electrode and PDMS lead colours.
for j=1:numel(keys) % Decode the same mesh surface as the browser.
    part=data.views{2}.(keys(j)); % Use station one in source order.
    vertices=reshape(typecast(matlab.net.base64decode(part.xyz),'single'),3,[])'; % Read world coordinates in mm.
    faces=double(reshape(typecast(matlab.net.base64decode(part.tri),'uint32'),3,[])')+1; % Convert triangle indices to MATLAB.
    patch(ax,'Vertices',vertices,'Faces',faces,'FaceColor',colours(j,:),'EdgeColor',[0.25 0.23 0.21],'EdgeAlpha',0.34,'LineWidth',0.35); % Draw the actual exposed tetrahedral faces.
end % Finish the selected meshes.
axis(ax,'equal','tight'); view(ax,[1 1.2 1.5]); xlabel(ax,'x (mm)');ylabel(ax,'y (mm)');zlabel(ax,'Along nerve (mm)'); % Keep equal spatial scales.
title(ax,sprintf('3D mesh at station %+d · epineurium cut away',sections(1))); % Identify the opening and its location.
legend(ax,{'Epineurium','Fascicles','Electrodes','PDMS lead'},'Location','eastoutside'); % Identify the visible materials.
exportgraphics(fig,fullfile(folder,'Cuff60_mesh_3d.png'),'Resolution',150);close(fig); % Save the gallery view for future reruns.
end
function patchData=surfaceData(X,T,omit) % Obtain triangular boundaries without smoothing or decimating.
[used,~,compact]=unique(T(:)); % Retain the referenced nodes of this material and window.
tri=used(freeBoundary(triangulation(reshape(compact,size(T)),X(used,:)))); % Remove faces shared by two selected tetrahedra.
if nargin>2 && ~isempty(omit), tri=tri(~ismember(sort(tri,2),omit,'rows'),:); end % Keep shared contact triangles on the tissue side only.
[nodes,~,idx]=unique(tri(:)); % Compact the original node indices.
patchData.xyz=matlab.net.base64encode(typecast(reshape(single(X(nodes,:)'),[],1),'uint8')); % Store physical xyz triples as float32.
patchData.tri=matlab.net.base64encode(typecast(reshape(uint32(reshape(idx,size(tri))'-1),[],1),'uint8')); % Encode zero-based vertex indices for WebGL.
end

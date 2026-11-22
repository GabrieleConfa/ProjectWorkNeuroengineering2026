function plotPersonalMesh(o,folder) % Draw a real section of the solved tetrahedra and their shared cuff face.
%% 1. Intersect the 3D tissue mesh with the first cuff station
s=o.sites(find(o.sites.type=="cuff",1),:); z=s.z_m+1e-9; % Offset the section by one nanometre to avoid coplanar vertices.
X=o.mesh.nodes_m; T=double(o.mesh.tets); cells=double(o.mesh.elemCell); tissue=o.cells.tissue(cells); tissue=tissue(:); % Access the mesh that produced the displayed field.
zz=reshape(X(3,T(:)),size(T)); cut=min(zz,[],1)<z & max(zz,[],1)>z; % Select each tetrahedron intersecting the transverse plane.
ids=T(:,cut); heights=zz(:,cut); n=size(ids,2); nearEdges=nchoosek(1:4,2); % Prepare the six possible segment intersections per tetrahedron.
px=nan(n,6); py=px; % Store a triangle or quadrilateral in each row.
for k=1:6 % Interpolate the intersection point on each crossing edge.
    a=nearEdges(k,1); b=nearEdges(k,2); hit=(heights(a,:)<z)~=(heights(b,:)<z); % Ignore non-crossing edges.
    t=(z-heights(a,hit))./(heights(b,hit)-heights(a,hit)); % Measure position along the edge.
    p=X(1:2,ids(a,hit))+(X(1:2,ids(b,hit))-X(1:2,ids(a,hit))).*t; % Calculate physical transverse coordinates.
    px(hit,k)=p(1,:)'*1e3; py(hit,k)=p(2,:)'*1e3; % Store millimetres for plotting.
end % Finish computing the cross-section polygons.
angle=atan2(py-mean(py,2,'omitnan'),px-mean(px,2,'omitnan')); [~,order]=sort(angle,2); % Sort valid vertices around each polygon centre.
idx=(1:n)'+(order-1)*n; px=px(idx); py=py(idx); % Match the connectivity and coordinates.
fig=figure(Visible="off",Theme="light",Color="white",Position=[80 80 1200 560]); layout=tiledlayout(fig,1,2); % Place the anatomical section next to the finite contact footprint.
ax=nexttile(layout);hold(ax,'on'); colours=[0.96 0.75 0.68;0.65 0.40 0.35;0.97 0.91 0.77;0.85 0.93 1;0.57 0.14 0.83]; % Use the tissue and electrode colours of the atlas.
for kind=unique(tissue(cut))' % Plot tissue groups with the same physical edge scale.
    for m=[3 4] % Each clipped tetrahedron has either three or four vertices.
        rows=tissue(cut)==kind & sum(isfinite(px),2)==m; % Select closed polygons of this kind.
        if any(rows),patch(ax,px(rows,1:m)',py(rows,1:m)',colours(kind,:),EdgeColor=[0.48 0.45 0.41],LineWidth=0.15);end % Draw true tetrahedral intersections.
    end % Finish both polygon sizes.
end % Finish all tissue regions.
axis(ax,'equal'); lim=1e3*o.patient.nerveRadius_m*1.12; xlim(ax,[-lim lim]);ylim(ax,[-lim lim]); % Frame the electrode station at its native aspect ratio.
xlabel(ax,'x (mm)');ylabel(ax,'y (mm)');title(ax,sprintf('Tetrahedra at station %+d · z = %.1f mm',s.section,z*1e3)); % Identify the slice.
%% 2. Display only triangles shared by the pad and the epineurium
[audit,faces]=checkCuffMeshContact(X,T,tissue,o.cells.site(cells),o.sites,o.params.personalGeometry); % Use the same interface test that gates the FEM solve.
j=find(audit.site==s.id,1); [nodes,~,index]=unique(faces{j}(:)); tri=reshape(index,size(faces{j})); % Compress the triangle indexing for this contact.
local=(X(:,nodes)'-s.centre_m)*1e3; uv=[local(:,1:2)*s.w(1:2)',local(:,3)]; % Rotate into the true transverse and axial pad coordinates.
ax=nexttile(layout);patch(ax,Faces=tri,Vertices=uv,FaceColor=[0.76 0.59 0.89],EdgeColor=[0.28 0.15 0.39]); % Show the actual connected interface, viewed normally.
axis(ax,'equal');grid(ax,'on');xlabel(ax,'Across pad (mm)');ylabel(ax,'Along nerve (mm)');title(ax,sprintf('AS%d · shared contact area %.4f mm²',s.id,audit.shared_mm2(j))); % Display physical area from the mesh.
title(layout,sprintf('3D conduction mesh · %s tetrahedra',string(size(T,2)))); % Report the full volume element count.
exportgraphics(fig,fullfile(folder,'Cuff60_mesh.png'),Resolution=155);close(fig); % Save the audited mesh for the personal gallery.
end

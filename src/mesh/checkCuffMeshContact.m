function [report,footprints]=checkCuffMeshContact(X,T,tissue,siteIds,sites,q) % Measure shared pad–epineurium triangles in the actual conduction mesh.
% Purpose: audit of the Gmsh mesh before solving. Each rigid cuff pad must share its whole active face with the epineurium
% (conforming triangles on the pad plane), otherwise the injected current would see a smaller contact than modelled. For each
% cuff pad it finds the triangles shared by the pad tetrahedra and the epineurium tetrahedra that lie on the pad plane
% (cuffContactPlane) and compares their area with the expected footprint (contact width x sideZ_m).
% Inputs:  X (3 x nNodes [m]), T (4 x nTets node indices), tissue (1 x nTets tissue code, 3 = epineurium), siteIds (1 x nTets
%          contact id of each element, 0 = none), sites (contact table), q (Gmsh-model parameters for cuffContactPlane).
% Outputs: report table per cuff pad (site, shared_mm2, expected_mm2, sharedTriangles, planeError_um); footprints (cell of
%          shared triangle node indices, for the figures). Errors if a pad is detached or its area differs by > 0.1 %.
cuff=sites(sites.type=="cuff",:); count=height(cuff); % Inspect every external electrode using its own local frame.
footprints=cell(count,1); % Retain the actual interface faces for the viewer and mesh figures.
area=zeros(count,1); expected=area; residual=area; triangles=area; % Store footprint and coplanarity diagnostics per electrode.
epi=T(:,tissue==3); % Restrict the interface search to epineurial tetrahedra.
for k=1:count % Audit each pad separately so an aggregate area cannot mask a gap.
    s=cuff(k,:); [plane,~,width,~,frame]=cuffContactPlane(q,s); % Read the prescribed inner plane and expected contact width.
    expected(k)=width*s.sideZ_m; n=frame.n'; % Define the expected active footprint in square metres.
    pad=T(:,siteIds==s.id); nodes=false(1,size(X,2)); nodes(pad(:))=true; % Identify the conductor nodes of this site.
    neighbours=epi(:,sum(nodes(epi),1)>=3); % Only an epineurial tetrahedron sharing three nodes can share a triangle.
    assert(~isempty(neighbours),'AS%d has no meshed epineurial interface.',s.id); % Reject a detached modelled electrode.
    padFaces=tetraFaces(pad); epiFaces=tetraFaces(neighbours); % Use global indices to match triangles across distinct tissue regions.
    shared=intersect(padFaces,epiFaces,'rows'); % Require identical triangles in both tissues.
    distance=abs(n'*X(1:2,:)-plane); inner=shared(all(distance(shared)<1e-10,2),:); % Select the rigid inner face separately from its side walls.
    assert(~isempty(inner),'AS%d has no shared triangle on the active face.',s.id); % Fail before solving a disconnected contact.
    a=X(:,inner(:,2))-X(:,inner(:,1)); b=X(:,inner(:,3))-X(:,inner(:,1)); % Form oriented triangle edges.
    footprints{k}=inner; % Store global triangle indices on this shared active face.
    area(k)=sum(vecnorm(cross(a,b,1),2,1))/2; residual(k)=max(distance(inner),[],'all'); triangles(k)=size(inner,1); % Integrate the matched footprint.
end % Finish all active faces.
report=table(cuff.id,area*1e6,expected*1e6,triangles,residual*1e6,'VariableNames',{'site','shared_mm2','expected_mm2','sharedTriangles','planeError_um'}); % Keep contact diagnostics with the numerical result.
assert(all(abs(area-expected)<1e-12+1e-3*expected),'Meshed cuff contact area differs from its geometric footprint.'); % Reject a partial interface in the default full-contact case.
end
function faces=tetraFaces(T) % Enumerate unoriented triangle keys for the local tetrahedra.
faces=sort([T([1 2 3],:) T([1 2 4],:) T([1 3 4],:) T([2 3 4],:)]',2); % Only keys present in both tissues can belong to the interface.
end

function X=fitCuffMeshSurface(X,T,cells,sites,geom,q) % Recover the local contacted shell from the landing-face reference extrusion.
%% 1. Intersect the actual node rays with the reference polygon
xy=geom.cuffReferenceXY; angles=mod(atan2(xy(:,2),xy(:,1)),2*pi); [angles,order]=sort(angles); xy=xy(order,:); % Sort irregular corners by polar angle.
band=0.12e-3; r=hypot(X(1,:),X(2,:)); dirs=X(1:2,:)'./max(r',eps); native=nerveBoundaryRadius(q,dirs); % Locate the original nerve boundary.
near=find(abs(r'-native)<band*1.1); fixed=false(1,size(X,2)); % Restrict motion to the exterior boundary layer.
fixed(T(:,cells==1 | cells==5))=true; near=near(~fixed(near)); % Preserve the fascicles and rigid electrode volumes exactly.
ray=dirs(near,:); a=mod(atan2(ray(:,2),ray(:,1))-angles(1),2*pi); % Locate one polygon edge for each candidate node.
j=discretize(a,[angles-angles(1);2*pi]); left=xy(j,:); right=xy(mod(j,size(xy,1))+1,:); edge=right-left; % Read the segment bracketing this node.
ref=(left(:,1).*edge(:,2)-left(:,2).*edge(:,1))./(ray(:,1).*edge(:,2)-ray(:,2).*edge(:,1)); % Intersect its ray with the actual reference edge.
%% 2. Shift the outer layer to the local pad footprint and smooth shoulders
cuff=sites(sites.type=="cuff",:); target=cuffDeformedRadius(q,cuff,X(3,near),ray,true)'; % Sample the same surface used by the viewer at each mesh node.
fraction=min(1,abs(r(near)'-ref)/band); weight=1-fraction.^2.*(3-2*fraction); % Extend the correction smoothly into the adjoining tissue.
original=X; % Keep the reference positions for signed volume checks.
shift=(target-ref).*weight; X(1:2,near)=X(1:2,near)+(shift.*ray)'; % Preserve axial coordinates and the physical pad faces.
%% 3. Reject any inverted or collapsed tetrahedra before assembling the field
changed=false(1,size(X,2)); changed(near)=abs(shift)>eps; % Locate only elements affected by the warp.
elements=find(any(changed(T),1)); % Unmoved tetrahedra retain their original volume.
for first=1:50000:numel(elements) % Limit transient memory on the Mac.
    ids=T(:,elements(first:min(first+49999,end))); % Gather a block of modified tetrahedra.
    ratio=signedVolume(X,ids)./signedVolume(original,ids); % Positive ratios preserve orientation.
    assert(all(isfinite(ratio) & ratio>0),'The contact fit inverted or collapsed a tetrahedron.'); % Reject inconsistent conduction geometry.
end % Finish checking the surface deformation.
end
function v=signedVolume(X,ids) % Return six times the oriented tetrahedral volume.
a=X(:,ids(1,:)); b=X(:,ids(2,:))-a; c=X(:,ids(3,:))-a; d=X(:,ids(4,:))-a; % Read the three edge vectors.
v=sum(cross(b,c,1).*d,1); % Preserve the orientation of each reference tetrahedron.
end

function [tetId,bary] = locateRegularGrid(X,T,axes_m,mask) % Locate regular voxels by testing tetrahedral bounding boxes.
%% 1. Prepare grid coordinates and output arrays
dims = cellfun(@numel,axes_m); % Count regular samples per spatial direction.
origin = cellfun(@(v)v(1),axes_m); % Read the first grid point in each direction.
step = cellfun(@(v)v(2)-v(1),axes_m); % Read uniform physical grid spacings.
ids = nan(prod(dims),1); weights = nan(prod(dims),4); % Keep points outside the mesh unassigned.
%% 2. Test only grid points inside each element's bounding box
for first=1:20000:size(T,2) % Bound temporary memory on the laptop.
    elements=first:min(first+19999,size(T,2)); % Select a contiguous set of tetrahedra.
    corners=T(:,elements); % Read the four vertex indices per tetrahedron.
    low=zeros(numel(elements),3); high=low; % Allocate inclusive bounding indices.
    for d=1:3 % Find voxel indices covered along each direction.
        coord=reshape(X(d,corners),4,[]); % Group the four vertex coordinates.
        low(:,d)=max(1,ceil((min(coord,[],1)'-origin(d))/step(d)-1e-10)+1); % Include lower shared faces.
        high(:,d)=min(dims(d),floor((max(coord,[],1)'-origin(d))/step(d)+1e-10)+1); % Include upper shared faces.
    end % Finish the three bounding intervals.
    span=max(0,high-low+1); counts=prod(span,2); % Count candidate voxels for each element.
    row=repelem((1:numel(elements))',counts); % Associate each candidate with its element.
    if isempty(row), continue; end % Skip batches outside the grid.
    row=row(:); starts=repelem(cumsum(counts)-counts,counts); offset=(0:numel(row)-1)'-starts(:); % Number candidates inside each bounding box.
    ijk=low(row,:)+[mod(offset,span(row,1)),mod(floor(offset./span(row,1)),span(row,2)),floor(offset./(span(row,1).*span(row,2)))]; % Decode the 3D indices.
    k=sub2ind(dims,ijk(:,1),ijk(:,2),ijk(:,3)); % Use MATLAB voxel order.
    keep=mask(k) & isnan(ids(k)); % Skip exterior and previously located points.
    k=k(keep); row=row(keep); ijk=ijk(keep,:); % Retain the necessary candidates only.
    pts=origin+(ijk-1).*step; % Convert candidate indices to physical coordinates.
    tc=corners(:,row); p1=X(:,tc(1,:)); % Locate the first vertex of each candidate element.
    a=X(:,tc(2,:))-p1; b=X(:,tc(3,:))-p1; c=X(:,tc(4,:))-p1; % Form the three edges.
    delta=pts'-p1; denom=dot(a,cross(b,c,1),1); % Form signed tetrahedral volumes.
    w=[dot(delta,cross(b,c,1),1);dot(delta,cross(c,a,1),1);dot(delta,cross(a,b,1),1)]./denom; % Compute three barycentric coordinates.
    w=[1-sum(w,1);w]'; % Recover the first coordinate by partition of unity.
    inside=all(w>=-1e-10 & w<=1+1e-10,2); % Accept points inside or on the tetrahedron.
    ids(k(inside))=elements(row(inside)); weights(k(inside),:)=w(inside,:); % Store matching element and linear interpolation weights.
end % Finish searching all physical elements.
tetId=ids(mask); bary=weights(mask,:); % Return the ordering expected by the field exporter.
end % Finish regular-grid location.

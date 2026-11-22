function xy=cuffReferenceSection(q,cuff) % Build a polygon with a full, exact landing face for each rigid pad orientation.
angle=linspace(0,2*pi,65)'; % Seed the rest of the circumference with regular vertices.
for k=1:height(cuff) % Include both finite edges of every landing face.
    s=cuff(k,:); [plane,~,~,~,frame]=cuffContactPlane(q,s); % Use the prescribed plane in the contact frame.
    angle=[angle;mod(frame.cornerAngles,2*pi)]; %#ok<AGROW> % Include the rays through both corners.
end % End edge sampling.
angle=uniquetol(angle,1e-12); angle(end)=[]; % Drop the duplicate closing direction.
ray=[cos(angle) sin(angle)]; radius=nerveBoundaryRadius(q,ray); % Locate the native ellipse on each ray.
for k=1:height(cuff) % Replace each short native arc with a planar patch for the reference extrusion.
    s=cuff(k,:); [plane,~,~,~,frame]=cuffContactPlane(q,s); % Read the physical normal and elevation.
    n=frame.n'; dotN=ray*n; t=ray*[-n(2);n(1)]; hit=dotN>0 & abs(plane*t./dotN-frame.offset)<=s.sideXY_m/2+1e-12; % Select only points inside the face.
    radius(hit)=min(radius(hit),plane./dotN(hit)); % Preserve the finite flat inner face.
end % Finish the reference section.
xy=radius.*ray; before=xy-circshift(xy,1); after=circshift(xy,-1)-xy; % Detect intermediate points on an exact straight chord.
keep=abs(before(:,1).*after(:,2)-before(:,2).*after(:,1))>1e-10*vecnorm(before,2,2).*vecnorm(after,2,2); % Retain only actual corners.
xy=xy(keep,:); % Return a counterclockwise, open list of polygon vertices.
end

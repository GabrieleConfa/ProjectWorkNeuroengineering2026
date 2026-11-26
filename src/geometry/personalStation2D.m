function g=personalStation2D(section,q) % Build the current cross-section for both of Gabriele's 2D workflows.
if nargin<2, q=quicklook3DParams(); end % Load current device and anatomy settings on a standalone run.
[~,~,sites]=personalDevice3D(pw6_3d_params(),q); % Recover the same finite pad poses as the interactive viewer.
g.sites=sites(sites.section==section,:); % Use the six physical contacts belonging to this section.
raw=loadActiveSites(q.electrode.sitesFile,q.geometry.radius_m); % Recover the incoming geometry for the approved lead planner.
cuff=raw(raw.section==section & raw.type=="cuff",:); % Select both opposing pads.
strip=planQuicklookStrip(q,raw(raw.section==section & raw.type=="intraneural",:),cuff,"deflect"); % Recover the continuous one-sided lead.
g.path=strip.path_m(:,1:2); g.f=q.geometry.fascicles_m; g.polys=stationFascicles(q,section); % Store the path, the legacy ellipses and the fascicle outlines of this station for the plots.
shell=buildQuicklookNerve(q,cuff); % Use the exact angular samples through the pad edges.
[~,j]=min(abs(shell.z_m(:,1)-cuff.z_m(1))); % Sample the centre plane of the pad.
xy=[shell.x_m(j,1:end-1)' shell.y_m(j,1:end-1)']; % Retain one vertex per circumferential ray.
prev=xy-circshift(xy,1); after=circshift(xy,-1)-xy; % Identify redundant vertices inside planar faces.
keep=abs(prev(:,1).*after(:,2)-prev(:,2).*after(:,1))>1e-16.*vecnorm(prev,2,2); % Keep the exact corners of each flat face.
g.outline=xy(keep,:); g.q=q; % Retain the native millimetre-independent geometry for future data.
end

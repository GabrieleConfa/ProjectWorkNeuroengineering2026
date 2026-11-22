classdef testCuffTangent < matlab.unittest.TestCase % Check the off-axis case that a circular contact convention misses.
    methods (Test)
        function rotatedOvalAndCircle(tc) % Verify the native tangent and the loaded finite face independently.
            q=quicklook3DParams(); raw=loadActiveSites(q.electrode.sitesFile,q.geometry.radius_m); % Read a real oblique contact.
            site=raw(raw.id==11,:); ray=[cosd(site.theta_deg) sind(site.theta_deg)]; % Preserve the supplied device azimuth.
            for b=[1.2e-3 2e-3] % Exercise an oval and its circular limit.
                q.geometry.nerveSemiAxes_m=[2e-3 b]; q.geometry.nerveRotation_deg=23; % Test a rotated anatomical frame.
                [~,depth,width,~,f]=cuffContactPlane(q,site); % Get the frame used by graphics and the mesher.
                anchor=ray*nerveBoundaryRadius(q,ray); % Locate the original polar landmark independently.
                plus=ray*[cos(1e-5) sin(1e-5);-sin(1e-5) cos(1e-5)]; % Rotate the ray slightly forward.
                minus=ray*[cos(1e-5) -sin(1e-5);sin(1e-5) cos(1e-5)]; % Rotate the ray slightly backward.
                tangent=plus*nerveBoundaryRadius(q,plus)-minus*nerveBoundaryRadius(q,minus); % Estimate the native tangent by centred differencing.
                tc.verifyLessThan(abs(dot(f.n,tangent/norm(tangent))),1e-8); % Require orthogonality at the supplied landmark.
                tc.verifyEqual(f.centre,anchor-depth*f.n,'AbsTol',1e-14); % Check translation along the physical normal.
                tc.verifyEqual(width,site.sideXY_m,'AbsTol',1e-12); % Verify both active edges are seated.
                target=f.centre+[-0.5;0;0.5]*width*f.t; dirs=target./vecnorm(target,2,2); % Sample the two real pad edges and centre.
                radii=cuffDeformedRadius(q,site,site.z_m,dirs); % Obtain the contacted tissue radius from the surface generator.
                tc.verifyEqual(radii'.*dirs,target,'AbsTol',1e-12); % Require tissue and electrode coincidence over the full face.
            end % Finish the rotated ellipse and circle checks.
        end
    end
end

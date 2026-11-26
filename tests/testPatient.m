classdef testPatient < matlab.unittest.TestCase % The patient layer: virtual and patient-specific data are interchangeable.
    methods (Test)
        function referencePatientIsValidWithoutFunctions(tc) % V001: valid; functions are unknown by default, so none are shown and no motor patches exist.
            q = quicklook3DParams(); pt = q.geometry.patient; tc.verifyEqual(pt.id,'V001'); validatePatient(pt,q); % Does not error.
            tc.verifyEqual(pt.functionSource,'unknown'); tc.verifyTrue(all(arrayfun(@(f)isempty(f.composition) && isempty(f.patches),pt.fascicles)));
        end
        function functionHypothesisGivesConsistentPatches(tc) % With an explicit hypothesis: compositions sum to 1, motor fraction about 13 percent, patch area = motor fraction.
            q = quicklook3DParams(struct('functionScenario',"random")); pt = q.geometry.patient; validatePatient(pt,q);
            C = cell2mat({pt.fascicles.composition}'); tc.verifyEqual(sum(C,2),ones(10,1),'AbsTol',1e-9); mf = [pt.fascicles.motorFraction]; tc.verifyGreaterThan(mean(mf),0.05); tc.verifyLessThan(mean(mf),0.25);
            for k = 1:10, a = sum(cellfun(@(Q)polyarea(Q(:,1),Q(:,2)),pt.fascicles(k).patches)); tc.verifyEqual(a/pt.fascicles(k).area_m2,mf(k),'RelTol',0.08); end
        end
        function fileRoundTrip(tc) % Saving and loading keeps outlines, compositions and patches.
            q = quicklook3DParams(); pt = q.geometry.patient; f = [tempname '.json']; cleanup = onCleanup(@()delete(f)); savePatient(pt,f); b = loadPatient(f);
            for k = 1:10, tc.verifyEqual(b.fascicles(k).outline,pt.fascicles(k).outline,'AbsTol',1e-12); tc.verifyEqual(numel(b.fascicles(k).composition),numel(pt.fascicles(k).composition)); tc.verifyEqual(numel(b.fascicles(k).patches),numel(pt.fascicles(k).patches)); end
        end
        function specificPatientRunsThroughTheModel(tc) % A patient given as a file (here the template) gets leads on every station.
            q = quicklook3DParams(struct('name',"TEMPLATE_measured")); sc = q.geometry.scenario; tc.verifyEqual(numel(sc.stations),8); tc.verifyEqual(q.geometry.patient.source,'measured');
            for s = sc.stations, tc.verifyEqual(size(s.padXY,1),4); tc.verifyFalse(any(isnan(s.padXY(:)))); end
            tc.verifyEqual(q.geometry.patient.functionSource,'unknown'); % One fascicle has no function: the model still works.
        end
        function patientFromOutlinesAndErrors(tc) % Helper for segmented sections, and a clear error for a wrong function id.
            th = linspace(0,2*pi,25)'; th(end) = []; nerve = 2e-3*[cos(th) 0.6*sin(th)]; f1 = 0.4e-3*[cos(th) sin(th)]-[0.8e-3 0];
            pt = patientFromOutlines("X",nerve,{f1},14e-6,"AIN"); validatePatient(pt,quicklook3DParams()); tc.verifyEqual(find(pt.fascicles(1).composition==1),4); % AIN is the 4th group.
            tc.verifyError(@()patientFromOutlines("X",nerve,{f1},14e-6,"NOPE"),?MException);
        end
        function functionHypothesesAreReproducible(tc) % Same seed, same assignment; different hypotheses differ.
            q = quicklook3DParams(); pt = q.geometry.patient; a = assignFunctions(pt,"random",5); b = assignFunctions(pt,"random",5); c = assignFunctions(pt,"mrAIN",5);
            tc.verifyEqual(cell2mat({a.fascicles.composition}'),cell2mat({b.fascicles.composition}')); tc.verifyNotEqual(cell2mat({a.fascicles.composition}'),cell2mat({c.fascicles.composition}'));
        end
    end
end

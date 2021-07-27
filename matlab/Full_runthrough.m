clear;

% geometrical parameters
x_gap = 20;                 % gap between emitter and collector in the x direction 
y_gap = 50;                % gap between 2 gates in the y direction
w = 100;                    % width of emitter and collector
wg = 80;                    % width of gates
l = 100;                     % length of emitter
lg = 100;                   % length of gate
h = 50;                     % height of triangle
hg = 50;                    % height of triangle of gate
r_e = 5;                    % radius of emmiter tip nm
r_g = 10;                   % radius of gate tip]
tol = 0.05;

Emitter = 'Emitter';
Collector = 'Collector';
Gate = 'Gate';

% set up general shape
emitter_vertices = [-x_gap/2, 0;-x_gap/2-h, w/2;-x_gap/2-h-l, w/2;-x_gap/2-h-l, 0;-x_gap/2, 0];
collector_vertices = [x_gap/2, 0; x_gap/2+h+l, 0; x_gap/2+h+l, w/2; x_gap/2+h, w/2;x_gap/2, 0];
gate_vertices = [0, y_gap/2;wg/2, y_gap/2+hg; wg/2,  y_gap/2+hg+lg; -wg/2,  y_gap/2+hg+lg; -wg/2, y_gap/2+hg; 0, y_gap/2;];

template = 'C:\Users\John\Desktop\Lorentz_Scripted\template.tpl';
% initialization
Lorentz = actxserver('IES.document');
iErr = 1;
while iErr > 0      % pause until Lorentz is completely initialized
  [iErr] = Window_CheckInitialization(Lorentz, iErr);  
end
[iErr] = File_OpenTemplate(Lorentz,template, iErr);     %load template
Model_Delete_All(Lorentz); %clear everything to begin


%set nm as the unit
[iErr] = Model_SetUnit(Lorentz, 'Length', 'nm', iErr);
[iErr] = View_Set2DLimits(Lorentz, -2*100, -2*100, 2*100, 2*100, iErr);

% set_up_geometry crates the segments for each element and rounds the corners 
[emitter_points, emitterSegs, collector_points, collectorSegs, gate_points, gateSegs] = set_up_geometry(Lorentz, emitter_vertices, collector_vertices, gate_vertices, Emitter, Collector, Gate, iErr);

% visualize geometry
figure(1)
plot(emitter_points(:,1),emitter_points(:,2), '.-',collector_points(:,1),collector_points(:,2), '.-', gate_points(:,1),gate_points(:,2), '.-')
% Enhanced Fowler-Nordheim parameters
Work_Fun = 4.6;
beta = 100;
alpha = 1e-8;

Ve = 0;
Vg = linspace(0,10,3);
Vc = linspace(0,15,3);
numPointsReq = 2;

gate_current = zeros(length(Vc), length(Vg));
collector_current = zeros(length(Vc), length(Vg));
for j=1:numel(Vg)
    for i=1:numel(Vc)
        [gate_current(i, j), collector_current(i, j)] = getCurrents_modified(Lorentz, Ve, Vc(i), Vg(j), Emitter, Gate, Collector, Work_Fun, beta, alpha, emitterSegs, gateSegs, collectorSegs, collector_points, gate_points, numPointsReq, 0, iErr);
    end
end


score = scoreGeom(Ve, Vg, Vc, collector_current, gate_current);
seg_len = 5;
NumPoints = length(emitter_points) + length(gate_points) + length(collector_points);
n = round(sqrt(NumPoints));
score_arr = score;
num_perts = 0;
iternums = 0;
for iterNum = 1:3
    inds = randi(NumPoints, 1, n);
    inds = unique(inds);
    for k= inds 
        num_perts = num_perts + 1;
        new_emitter_points = emitter_points;
        new_collector_points = collector_points;
        new_gate_points = gate_points;
        if k < length(emitter_points)+1
            % change k and k +1 in emitter
            q = k;
            pert = getPerturbation(emitter_points, q, 2*seg_len, tol);
            if q == 1 || q == length(emitter_points)
                new_emitter_points(1, :) = new_emitter_points(1, :)+ pert;
                new_emitter_points(end, :) = new_emitter_points(end, :)+ pert;
                [iErr] = Object_RemoveSegment(Lorentz, Emitter, emitterSegs(1), iErr);
                [iErr] = Object_RemoveSegment(Lorentz, Emitter, emitterSegs(end), iErr);
                % apply pert to emitter
                [emitterSegs(end), iErr] = Geometry2D_CreateLine(Lorentz, new_emitter_points(end-1,1), new_emitter_points(end-1,2), new_emitter_points(end,1), new_emitter_points(end,2), emitterSegs(end), iErr);
                [emitterSegs(1), iErr] = Geometry2D_CreateLine(Lorentz, new_emitter_points(1,1), new_emitter_points(1,2),new_emitter_points(2,1), new_emitter_points(2,2), emitterSegs(1), iErr);
                [iErr] = Object_AddSegment(Lorentz, Emitter, emitterSegs(end), iErr);
                [iErr] = Object_AddSegment(Lorentz, Emitter, emitterSegs(1), iErr);
                % evaluate effect of change 
                for j=1:numel(Vg)
                    for i=1:numel(Vc)
                        [gate_current(i, j), collector_current(i, j)] = getCurrents_modified(Lorentz, Ve, Vc(i), Vg(j), Emitter, Gate, Collector, Work_Fun, beta, alpha, emitterSegs, gateSegs, collectorSegs, new_collector_points, new_gate_points, numPointsReq, num_perts, iErr);
                    end
                end
                new_score = scoreGeom(Ve, Vg, Vc, collector_current, gate_current);
                if new_score < score
                    % keep the changes
                    score = new_score;
                    emitter_points = new_emitter_points;
                else
                    % discard the changes
                    [iErr] = Object_RemoveSegment(Lorentz, Emitter, emitterSegs(1), iErr);
                    [iErr] = Object_RemoveSegment(Lorentz, Emitter, emitterSegs(end), iErr);
                    [emitterSegs(end), iErr] = Geometry2D_CreateLine(Lorentz, emitter_points(end-1,1), emitter_points(end-1,2), emitter_points(end,1), emitter_points(end,2), emitterSegs(end), iErr);
                    [emitterSegs(1), iErr] = Geometry2D_CreateLine(Lorentz, emitter_points(1,1), emitter_points(1,2), emitter_points(2,1), emitter_points(2,2), emitterSegs(1), iErr);
                    [iErr] = Object_AddSegment(Lorentz, Emitter, emitterSegs(end), iErr);
                    [iErr] = Object_AddSegment(Lorentz, Emitter, emitterSegs(1), iErr);
                end
            else
                new_emitter_points(q, :) = new_emitter_points(q, :)+ pert;
                [iErr] = Object_RemoveSegment(Lorentz, Emitter, emitterSegs(q-1), iErr);
                [iErr] = Object_RemoveSegment(Lorentz, Emitter, emitterSegs(q), iErr);
                % apply pert to emitter
                [emitterSegs(q-1), iErr] = Geometry2D_CreateLine(Lorentz, new_emitter_points(q -1,1), new_emitter_points(q-1,2), new_emitter_points(q,1), new_emitter_points(q,2), emitterSegs(q-1), iErr);
                [emitterSegs(q), iErr] = Geometry2D_CreateLine(Lorentz, new_emitter_points(q,1), new_emitter_points(q,2),new_emitter_points(q+1,1), new_emitter_points(q+1,2), emitterSegs(q), iErr);
                [iErr] = Object_AddSegment(Lorentz, Emitter, emitterSegs(q-1), iErr);
                [iErr] = Object_AddSegment(Lorentz, Emitter, emitterSegs(q), iErr);
                % evaluate effect of change
                for j=1:numel(Vg)
                    for i=1:numel(Vc)
                        [gate_current(i, j), collector_current(i, j)] = getCurrents_modified(Lorentz, Ve, Vc(i), Vg(j), Emitter, Gate, Collector, Work_Fun, beta, alpha, emitterSegs, gateSegs, collectorSegs, new_collector_points, new_gate_points, numPointsReq, num_perts, iErr);
                    end
                end
                new_score = scoreGeom(Ve, Vg, Vc, collector_current, gate_current);
                if new_score < score
                    % keep the changes
                    score = new_score;
                    emitter_points = new_emitter_points;
                else
                    % discard the changes
                    [iErr] = Object_RemoveSegment(Lorentz, Emitter, emitterSegs(q-1), iErr);
                    [iErr] = Object_RemoveSegment(Lorentz, Emitter, emitterSegs(q), iErr);
                    [emitterSegs(q-1), iErr] = Geometry2D_CreateLine(Lorentz, emitter_points(q -1,1), emitter_points(q-1,2), emitter_points(q,1), emitter_points(q,2), emitterSegs(q-1), iErr);
                    [emitterSegs(q), iErr] = Geometry2D_CreateLine(Lorentz, emitter_points(q,1), emitter_points(q,2), emitter_points(q+1,1), emitter_points(q+1,2), emitterSegs(q), iErr);
                    [iErr] = Object_AddSegment(Lorentz, Emitter, emitterSegs(q-1), iErr);
                    [iErr] = Object_AddSegment(Lorentz, Emitter, emitterSegs(q), iErr);
                end
            end
        elseif k < length(collector_points) + length(emitter_points) +1
            % change j - emitter.length and j - emitter.length +1 in
            % collector
            q = k - length(emitter_points);
            pert = getPerturbation(collector_points, q, 2*seg_len, tol);
            if q == 1 || q == length(collector_points)
                new_collector_points(1, :) = new_collector_points(1, :) + pert;
                new_collector_points(end, :) = new_collector_points(end, :) + pert;
                [iErr] = Object_RemoveSegment(Lorentz, Collector, collectorSegs(end), iErr);
                [iErr] = Object_RemoveSegment(Lorentz, Collector, collectorSegs(1), iErr);
                % apply pert to collector
                [collectorSegs(end), iErr] = Geometry2D_CreateLine(Lorentz, new_collector_points(end-1,1), new_collector_points(end-1,2), new_collector_points(end,1), new_collector_points(end,2), collectorSegs(end), iErr);
                [collectorSegs(1), iErr] = Geometry2D_CreateLine(Lorentz, new_collector_points(1,1), new_collector_points(1,2),new_collector_points(2,1), new_collector_points(2,2), collectorSegs(1), iErr);
                [iErr] = Object_AddSegment(Lorentz, Collector, collectorSegs(end), iErr);
                [iErr] = Object_AddSegment(Lorentz, Collector, collectorSegs(1), iErr);
                % evaluate effect of change
                for j=1:numel(Vg)
                    for i=1:numel(Vc)
                        [gate_current(i, j), collector_current(i, j)] = getCurrents_modified(Lorentz, Ve, Vc(i), Vg(j), Emitter, Gate, Collector, Work_Fun, beta, alpha, emitterSegs, gateSegs, collectorSegs, new_collector_points, new_gate_points, numPointsReq, num_perts, iErr);
                    end
                end
                new_score = scoreGeom(Ve, Vg, Vc, collector_current, gate_current);
                if new_score < score
                    % keep changes
                    score = new_score;
                    collector_points = new_collector_points;
                else
                    % discard changes
                    [iErr] = Object_RemoveSegment(Lorentz, Collector, collectorSegs(end), iErr);
                    [iErr] = Object_RemoveSegment(Lorentz, Collector, collectorSegs(1), iErr);
                    [collectorSegs(end), iErr] = Geometry2D_CreateLine(Lorentz, collector_points(end-1,1), collector_points(end-1,2), collector_points(end,1), collector_points(end,2), collectorSegs(end), iErr);
                    [collectorSegs(1), iErr] = Geometry2D_CreateLine(Lorentz, collector_points(1,1), collector_points(1,2), collector_points(2,1), collector_points(2,2), collectorSegs(1), iErr);
                    [iErr] = Object_AddSegment(Lorentz, Collector, collectorSegs(end), iErr);
                    [iErr] = Object_AddSegment(Lorentz, Collector, collectorSegs(1), iErr);
                end
            else
                new_collector_points(q, :) = new_collector_points(q, :)+ pert;
                [iErr] = Object_RemoveSegment(Lorentz, Collector, collectorSegs(q-1), iErr);
                [iErr] = Object_RemoveSegment(Lorentz, Collector, collectorSegs(q), iErr);
                % apply pert to collector
                [collectorSegs(q-1), iErr] = Geometry2D_CreateLine(Lorentz, new_collector_points(q -1,1), new_collector_points(q-1,2), new_collector_points(q,1), new_collector_points(q,2), collectorSegs(q-1), iErr);
                [collectorSegs(q), iErr] = Geometry2D_CreateLine(Lorentz, new_collector_points(q,1), new_collector_points(q,2),new_collector_points(q+1,1), new_collector_points(q+1,2), collectorSegs(q), iErr);
                [iErr] = Object_AddSegment(Lorentz, Collector, collectorSegs(q-1), iErr);
                [iErr] = Object_AddSegment(Lorentz, Collector, collectorSegs(q), iErr);
                % evaluate effect of change
                for j=1:numel(Vg)
                    for i=1:numel(Vc)
                        [gate_current(i, j), collector_current(i, j)] = getCurrents_modified(Lorentz, Ve, Vc, Vg, Emitter, Gate, Collector, Work_Fun, beta, alpha, emitterSegs, gateSegs, collectorSegs, new_collector_points, new_gate_points, numPointsReq, num_perts, iErr);
                    end
                end
                new_score = scoreGeom(Ve, Vg, Vc, collector_current, gate_current);
                if new_score < score
                    % keep the changes
                    score = new_score;
                    collector_points = new_collector_points;
                else
                    % discard the changes
                    [iErr] = Object_RemoveSegment(Lorentz, Collector, collectorSegs(q-1), iErr);
                    [iErr] = Object_RemoveSegment(Lorentz, Collector, collectorSegs(q), iErr);
                    % apply pert to collector
                    [collectorSegs(q-1), iErr] = Geometry2D_CreateLine(Lorentz, collector_points(q -1,1), collector_points(q-1,2), collector_points(q,1), collector_points(q,2), collectorSegs(q-1), iErr);
                    [collectorSegs(q), iErr] = Geometry2D_CreateLine(Lorentz, collector_points(q,1), collector_points(q,2),new_collector_points(q+1,1), collector_points(q+1,2), collectorSegs(q), iErr);
                    [iErr] = Object_AddSegment(Lorentz, Collector, collectorSegs(q-1), iErr);
                    [iErr] = Object_AddSegment(Lorentz, Collector, collectorSegs(q), iErr);
                end
            end
        else
            % change j - emitter.length - collector.lrngth in gate
            q = k - length(emitter_points) - length(collector_points);
            pert = getPerturbation(gate_points, q, 2*seg_len, tol);
            if q == 1 || q == length(gate_points)
                new_gate_points(1, :) = new_gate_points(1, :) + pert;
                new_gate_points(end, :) = new_gate_points(end, :) + pert;
                [iErr] = Object_RemoveSegment(Lorentz, Gate, gateSegs(end), iErr);
                [iErr] = Object_RemoveSegment(Lorentz, Gate, gateSegs(1), iErr);
                % apply pert to gate
                [gateSegs(end), iErr] = Geometry2D_CreateLine(Lorentz, new_gate_points(end-1,1), new_gate_points(end-1,2), new_gate_points(end,1), new_gate_points(end,2), gateSegs(end), iErr);
                [gateSegs(1), iErr] = Geometry2D_CreateLine(Lorentz, new_gate_points(1,1), new_gate_points(1,2),new_gate_points(2,1), new_gate_points(2,2), gateSegs(1), iErr);
                [iErr] = Object_AddSegment(Lorentz, Gate, gateSegs(end), iErr);
                [iErr] = Object_AddSegment(Lorentz, Gate, gateSegs(1), iErr);
                % evaluate effect of change
                for j=1:numel(Vg)
                    for i=1:numel(Vc)
                        [gate_current(i, j), collector_current(i, j)] = getCurrents_modified(Lorentz, Ve, Vc(i), Vg(j), Emitter, Gate, Collector, Work_Fun, beta, alpha, emitterSegs, gateSegs, collectorSegs, new_collector_points, new_gate_points, numPointsReq, num_perts, iErr);
                    end
                end
                new_score = scoreGeom(Ve, Vg, Vc, collector_current, gate_current);
                if new_score < score
                    score = new_score;
                    gate_points = new_gate_points;
                else
                    [iErr] = Object_RemoveSegment(Lorentz, Gate, gateSegs(end), iErr);
                    [iErr] = Object_RemoveSegment(Lorentz, Gate, gateSegs(1), iErr);
                    % apply pert to gate
                    [gateSegs(end), iErr] = Geometry2D_CreateLine(Lorentz, gate_points(end-1,1), gate_points(end-1,2), gate_points(end,1), gate_points(end,2), gateSegs(end), iErr);
                    [gateSegs(1), iErr] = Geometry2D_CreateLine(Lorentz, gate_points(1,1), gate_points(1,2), gate_points(2,1), gate_points(2,2), gateSegs(1), iErr);
                    [iErr] = Object_AddSegment(Lorentz, Gate, gateSegs(end), iErr);
                    [iErr] = Object_AddSegment(Lorentz, Gate, gateSegs(1), iErr);
                end   
            else
                new_gate_points(q, :) = new_gate_points(q, :)+ pert;
                [iErr] = Object_RemoveSegment(Lorentz, Gate, gateSegs(q-1), iErr);
                [iErr] = Object_RemoveSegment(Lorentz, Gate, gateSegs(q), iErr);
                % apply pert to gate
                [gateSegs(q-1), iErr] = Geometry2D_CreateLine(Lorentz, new_gate_points(q-1,1), new_gate_points(q-1,2), new_gate_points(q,1), new_gate_points(q,2), gateSegs(q-1), iErr);
                [gateSegs(q), iErr] = Geometry2D_CreateLine(Lorentz, new_gate_points(q,1), new_gate_points(q,2),new_gate_points(q+1,1), new_gate_points(q+1,2), gateSegs(q), iErr);
                [iErr] = Object_AddSegment(Lorentz, Gate, gateSegs(q-1), iErr);
                [iErr] = Object_AddSegment(Lorentz, Gate, gateSegs(q), iErr);
                % evaluate effect of change
                for j=1:numel(Vg)
                    for i=1:numel(Vc)
                        [gate_current(i, j), collector_current(i, j)] = getCurrents_modified(Lorentz, Ve, Vc(i), Vg(j), Emitter, Gate, Collector, Work_Fun, beta, alpha, emitterSegs, gateSegs, collectorSegs, new_collector_points, new_gate_points, numPointsReq, num_perts, iErr);
                    end
                end
                new_score = scoreGeom(Ve, Vg, Vc, collector_current, gate_current);
                if new_score < score
                    score = new_score;
                    gate_points = new_gate_points;
                else
                    [iErr] = Object_RemoveSegment(Lorentz, Gate, gateSegs(q-1), iErr);
                    [iErr] = Object_RemoveSegment(Lorentz, Gate, gateSegs(q), iErr);
                    % apply pert to gate
                    [gateSegs(q-1), iErr] = Geometry2D_CreateLine(Lorentz, gate_points(q-1,1), gate_points(q-1,2), gate_points(q,1), gate_points(q,2), gateSegs(q-1), iErr);
                    [gateSegs(q), iErr] = Geometry2D_CreateLine(Lorentz, gate_points(q,1), gate_points(q,2), gate_points(q+1,1), gate_points(q+1,2), gateSegs(q), iErr);
                    [iErr] = Object_AddSegment(Lorentz, Gate, gateSegs(q-1), iErr);
                    [iErr] = Object_AddSegment(Lorentz, Gate, gateSegs(q), iErr);
                end
            end
        end
        figure(1)
        plot(emitter_points(:,1),emitter_points(:,2), '.-',collector_points(:,1),collector_points(:,2), '.-', gate_points(:,1),gate_points(:,2), '.-');
        save(['Data/Geometry_iter_' num2str(num_perts) '.mat'], 'emitter_points', 'collector_points', 'gate_points')
        score_arr = [score_arr; score];
        iternums = [iternums; num_perts];
        figure(2)
        plot(iternums, score_arr)
    end
end
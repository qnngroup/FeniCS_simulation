function [emitter_points, emitter_seg_array, collector_points, collector_seg_array, gate_points, gate_seg_array] = set_up_geometry(Lorentz, emitter_vertices, collector_vertices, gate_vertices, Emitter, Collector, Gate, iErr)

Obj1 = int32(1);
Obj2 = int32(1);
Obj3 = int32(1);

 
seg_len = 7;

% emitter initial geometry
emitter_points = [;];
em_indecies_of_corners = [];
corner_ind = 1;
for i=1:length(emitter_vertices)-1
    x1 = emitter_vertices(i, 1); y1 = emitter_vertices(i, 2);
    x2 = emitter_vertices(i+1, 1); y2 = emitter_vertices(i+1, 2);
    points = getSegmentationBetweenPoints(x1, y1, x2, y2, seg_len);
    emitter_points = [emitter_points; points];
    em_indecies_of_corners = [em_indecies_of_corners; corner_ind];
    corner_ind = corner_ind + length(points);
end
emitter_points = [emitter_points; emitter_vertices(end, :)];
em_indecies_of_corners = [em_indecies_of_corners; corner_ind];

% round emitter corners
%figure;
for j = 2:length(em_indecies_of_corners)-1
    ind1 = em_indecies_of_corners(j) -3;
    ind2 = em_indecies_of_corners(j) +3;
    [pts, xs, ys] = roundCorners(emitter_points, ind1, ind2, 5, 0.01);
    %plot(xs, ys)
    %hold on
    emitter_points(ind1+1:ind2-1, :) = pts;
end
[pts, xs, ys] = roundCorners(emitter_points, length(emitter_points)-2, 3, 3, 0.01);
emitter_points_pr = emitter_points(3:end-2,:);
emitter_points = [emitter_points(end-2,:); pts; emitter_points_pr];
ind = find(emitter_points(:,2) ==0);
emitter_points = [emitter_points(1:ind(2), :); emitter_points(end, :)];
%plot(xs, ys)
%plot(emitter_vertices(:,1),emitter_vertices(:,2), 'o', emitter_points(:,1),emitter_points(:,2), '.-')


% Collector initial geometry 
collector_points = [;];
col_indecies_of_corners = [];
corner_ind = 1;
for i=1:length(collector_vertices)-1
    x1 = collector_vertices(i, 1); y1 = collector_vertices(i, 2);
    x2 = collector_vertices(i+1, 1); y2 = collector_vertices(i+1, 2);
    points = getSegmentationBetweenPoints(x1, y1, x2, y2, seg_len);
    collector_points = [collector_points; points];
    col_indecies_of_corners = [col_indecies_of_corners; corner_ind];
    corner_ind = corner_ind + length(points);
end
collector_points = [collector_points; collector_vertices(end, :)];
col_indecies_of_corners = [col_indecies_of_corners; corner_ind];

% Round collector corners
for j = 2:length(col_indecies_of_corners)-1
    ind1 = col_indecies_of_corners(j) -3;
    ind2 = col_indecies_of_corners(j) +3;
    [pts, xs, ys] = roundCorners(collector_points, ind1, ind2, 5, 0.01);
    %plot(xs, ys)
    %hold on
    collector_points(ind1+1:ind2-1, :) = pts;
end
[pts, xs, ys] = roundCorners(collector_points, length(collector_points)-2,3, 3, 0.01);
collector_points_pr = collector_points(3:end-2,:);
collector_points = [collector_points_pr; pts; collector_points(3, :)];
ind = find(collector_points(:,2) ==0);
collector_points = [collector_points(1, :); collector_points(ind(end-1):end, :)];
%plot(xs, ys)
%plot(collector_vertices(:,1),collector_vertices(:,2), 'o', collector_points(:,1),collector_points(:,2), '.-')

% gate initial geometry
gate_points = [;];
gt_indecies_of_corners = [];
corner_ind = 1;
for i=1:length(gate_vertices)-1
    x1 = gate_vertices(i, 1); y1 = gate_vertices(i, 2);
    x2 = gate_vertices(i+1, 1); y2 = gate_vertices(i+1, 2);
    points = getSegmentationBetweenPoints(x1, y1, x2, y2, seg_len);
    gate_points = [gate_points; points];
    gt_indecies_of_corners = [gt_indecies_of_corners; corner_ind];
    corner_ind = corner_ind + length(points);
end
gate_points = [gate_points;  gate_vertices(end, :)];
gt_indecies_of_corners = [gt_indecies_of_corners; corner_ind];

% Round collector corners
for j = 2:length(gt_indecies_of_corners)-1
    ind1 = gt_indecies_of_corners(j) -3;
    ind2 = gt_indecies_of_corners(j) +3;
    [pts, xs, ys] = roundCorners(gate_points, ind1, ind2, 5, 0.01);
    %plot(xs, ys)
    gate_points(ind1+1:ind2-1, :) = pts;
end
[pts, xs, ys] = roundCorners(gate_points, length(gate_points)-2,3, 5, 0.01);
gate_points_pr = gate_points(3:end-2,:);
gate_points = [gate_points_pr; pts; gate_points(3, :)]; 
%plot(gate_vertices(:,1), gate_vertices(:,2), 'o', gate_points(:,1), gate_points(:,2), '.-')

% Create emitter in lorentz
emitter_seg_array = [];
for i=1:length(emitter_points)-1
    seg = int32(1);
    [seg, iErr] = Geometry2D_CreateLine(Lorentz, emitter_points(i,1), emitter_points(i,2),emitter_points(i+1,1), emitter_points(i+1,2), seg, iErr);
    emitter_seg_array = [emitter_seg_array; seg];
end
% add emmiter segments to object
[Obj1, iErr] = Object_Create(Lorentz, Emitter, Obj1, iErr);
for i=1:length(emitter_seg_array)-1
    [iErr] = Object_AddSegment(Lorentz, Emitter, emitter_seg_array(i), iErr);
end 

% Create collector in lorentz
collector_seg_array = [];
for i=1:length(collector_points)-1
    seg = int32(1);
    [seg, iErr] = Geometry2D_CreateLine(Lorentz, collector_points(i,1), collector_points(i,2),collector_points(i+1,1), collector_points(i+1,2), seg, iErr);
    collector_seg_array = [collector_seg_array; seg];
end

% add collector segments to object
[Obj2, iErr] = Object_Create(Lorentz, Collector, Obj2, iErr);
for i=2:length(collector_seg_array)
    [iErr] = Object_AddSegment(Lorentz, Collector, collector_seg_array(i), iErr);   
end 

% Create gate in lorentz
gate_seg_array = [];
for i=1:length(gate_points)-1
    seg = int32(1);
    [seg, iErr] = Geometry2D_CreateLine(Lorentz, gate_points(i,1), gate_points(i,2), gate_points(i+1,1), gate_points(i+1,2), seg, iErr);
    gate_seg_array = [gate_seg_array; seg];
end 
% add gate segments to object
[Obj3, iErr] = Object_Create(Lorentz, Gate, Obj3, iErr);
for i=1:length(gate_seg_array)
    [iErr] = Object_AddSegment(Lorentz, Gate, gate_seg_array(i), iErr);   
end %




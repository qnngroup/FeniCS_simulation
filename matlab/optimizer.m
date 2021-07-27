clear;
close all;
%template with the right precision setup
template = 'C:\Users\John\Desktop\Lorentz_Scripted\template.tpl';

% Enhanced Fowler-Nordheim parameters
Work_Fun = 4.6;
beta = 100;
alpha = 1e-8;


% geometrical parameters
x_gap = 20;                 % gap between emitter and collector in the x direction 
y_gap = 50;                % gap between 2 gates in the y direction
w = 100;                    % width of emitter and collector
wg = 80;                    % width of gates
l = 200;                     % length of emitter
lg = 200;                   % length of gate
h = 50;                     % height of triangle
hg = 50;                    % height of triangle of gate
r_e = 5;                    % radius of emmiter tip nm
r_g = 10;                   % radius of gate tip
offsret = 0;      % offset in order to move things above the x axsis
r2d =180/pi;
Emitter = 'Emitter';
Collector = 'Collector';
Gate = 'Gate';
Obj1 = int32(1);
Obj2 = int32(1);
Obj3 = int32(1);


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
[iErr] = View_Set2DLimits(Lorentz, -3*w, -3*w, 3*w, 3*w, iErr);


% emitter initial geometry
seg_len = 5;

emitter_vertices = [-x_gap/2, 0;-x_gap/2-h, w/2;-x_gap/2-h-l, w/2;-x_gap/2-h-l, 0;-x_gap/2, 0];
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
emitter_points = [emitter_points; -x_gap/2, 0];

for j = 2:length(em_indecies_of_corners)-1
    ind1 = em_indecies_of_corners(j) -3;
    ind2 = em_indecies_of_corners(j) +3;
    pts = roundCorners(emitter_points, ind1, ind2, 5, 0.1);
    emitter_points(ind1+1:ind2-1, :) = pts;
end
pts = roundCorners(emitter_points, 2, length(emitter_points)-2, 3, 0.1);
emitter_points_pr = emitter_points(2:end-2,:);
emitter_points = [emitter_points(end-2,:); pts; emitter_points_pr];
ind = find(emitter_points(:,2) ==0);
emitter_points = [emitter_points(1:ind(2), :); emitter_points(end, :)];

num_segs = 0;
emitter_seg_array = [];
for i=1:length(emitter_points)-1
   
    seg = int32(1);
    [seg, iErr] = Geometry2D_CreateLine(Lorentz, emitter_points(i,1), emitter_points(i,2),emitter_points(i+1,1), emitter_points(i+1,2), seg, iErr);
    emitter_seg_array = [emitter_seg_array; seg];
    num_segs = num_segs+1;
end
seg = int32(1);
[seg, iErr] = Geometry2D_CreateLine(Lorentz, emitter_points(end,1), emitter_points(end,2),emitter_points(1,1), emitter_points(1,2), seg, iErr);
emitter_seg_array = [emitter_seg_array; seg];
%em_num_segs_array = [em_num_segs_array; 1];

[Obj1, iErr] = Object_Create(Lorentz, Emitter, Obj1, iErr);
for i=1:length(emitter_seg_array)-1
    [iErr] = Object_AddSegment(Lorentz, Emitter, emitter_seg_array(i), iErr);
end 

% collector initial geometry
collector_vertices = [x_gap/2, 0; x_gap/2+h+l, 0; x_gap/2+h+l, w/2; x_gap/2+h, w/2;x_gap/2, 0];
collector_points = [;];
col_num_points_on_seg = [];
for i=1:length(collector_vertices)-2
    x1 = collector_vertices(i, 1); y1 = collector_vertices(i, 2);
    x2 = collector_vertices(i+1, 1); y2 = collector_vertices(i+1, 2);
    points = getSegmentationBetweenPoints(x1, y1, x2, y2, seg_len);
    collector_points = [collector_points; points];
    col_num_points_on_seg = [col_num_points_on_seg; length(points)];
end

num_segs = 0;
col_num_segs_array = [];
collector_seg_array = [];
for i=1:length(collector_points)-1
    if collector_points(i, 1) ~= collector_points(i+1, 1) || collector_points(i, 2) ~= collector_points(i+1, 2)
        seg = int32(1);
        [seg, iErr] = Geometry2D_CreateLine(Lorentz, collector_points(i,1), collector_points(i,2),collector_points(i+1,1), collector_points(i+1,2), seg, iErr);
        collector_seg_array = [collector_seg_array; seg];
        num_segs = num_segs+1;
    else
        col_num_segs_array = [col_num_segs_array; num_segs];
        num_segs = 0;
    end 
end
seg = int32(1);
[seg, iErr] = Geometry2D_CreateLine(Lorentz, collector_points(end,1), collector_points(end,2),collector_points(1,1), collector_points(1,2), seg, iErr);
collector_seg_array = [collector_seg_array; seg];
col_num_segs_array = [col_num_segs_array; 1];


[Obj2, iErr] = Object_Create(Lorentz, Collector, Obj2, iErr);
for i=1:length(collector_seg_array)-1
    [iErr] = Object_AddSegment(Lorentz, Collector, collector_seg_array(i), iErr);   
end 


% gate initial geometry
gate_vertices = [0, y_gap/2;wg/2, y_gap/2+hg; wg/2,  y_gap/2+hg+lg; -wg/2,  y_gap/2+hg+lg; -wg/2, y_gap/2+hg; 0, y_gap/2;];
gate_points = [;];
gt_num_points_on_seg = [];
for i=1:length(gate_vertices)-1
    x1 = gate_vertices(i, 1); y1 = gate_vertices(i, 2);
    x2 = gate_vertices(i+1, 1); y2 = gate_vertices(i+1, 2);
    points = getSegmentationBetweenPoints(x1, y1, x2, y2, seg_len);
    gate_points = [gate_points; points];
    gt_num_points_on_seg = [gt_num_points_on_seg; length(points)];
end

num_segs = 0;
gt_num_segs_array = [];
gate_seg_array = [];
for i=1:length(gate_points)-1
    if gate_points(i, 1) ~= gate_points(i+1, 1) || gate_points(i, 2) ~= gate_points(i+1, 2)
        seg = int32(1);
        [seg, iErr] = Geometry2D_CreateLine(Lorentz, gate_points(i,1), gate_points(i,2), gate_points(i+1,1), gate_points(i+1,2), seg, iErr);
        gate_seg_array = [gate_seg_array; seg];
        num_segs = num_segs+1;
    else
        gt_num_segs_array = [gt_num_segs_array; num_segs];
        num_segs = 0;
    end 
end
gt_num_segs_array = [gt_num_segs_array; 1];

[Obj3, iErr] = Object_Create(Lorentz, Gate, Obj3, iErr);
for i=1:length(gate_seg_array)
    [iErr] = Object_AddSegment(Lorentz, Gate, gate_seg_array(i), iErr);   
end %


figure;
plot(emitter_vertices(:,1),emitter_vertices(:,2), 'o', emitter_points(:,1),emitter_points(:,2), '.-')
hold on 
plot(collector_vertices(:,1),collector_vertices(:,2), 'o', collector_points(:,1),collector_points(:,2), '.-')
plot(gate_vertices(:,1),gate_vertices(:,2), 'o', gate_points(:,1),gate_points(:,2), '.-')

% assign voltages 
% solve PDE
% collect data current 

Ve = 0;
Vg = linspace(0,10,4);
Vc = linspace(0,15,6);
numPointsReq = 3;

for j=1:numel(Vg)
    for i=1:numel(Vc)
    [gate_current(i, j), collector_current(i, j)] = getCurrents_modified(Lorentz, Ve, Vc(i), Vg(j), Emitter, Gate, Collector, Work_Fun, beta, alpha, emitter_seg_array, gate_seg_array, collector_seg_array, collector_points, gate_points, numPointsReq, iErr);
    end 
end

figure; 
plot(Vc,collector_current)
xlabel('Collector voltage')
ylabel('Collector current')

figure; 
plot(Vg,collector_current')
xlabel('Gate voltage')
ylabel('Collector current')

figure; 
plot(Vg,gate_current')
xlabel('gate voltage')
ylabel('gate current')

Window_Close(Lorentz)
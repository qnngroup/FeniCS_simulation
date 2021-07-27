function [emitter_seg_array, collector_seg_array, gate_seg_array] = re_set_up_geometry(Lorentz, emitter_points, collector_points, gate_points, Emitter, Collector, Gate, iErr)

Geometry_Delete_All(Lorentz);

% Create emitter in lorentz
emitter_seg_array = [];
for i=1:length(emitter_points)-1
    seg = int32(1);
    [seg, iErr] = Geometry2D_CreateLine(Lorentz, emitter_points(i,1), emitter_points(i,2),emitter_points(i+1,1), emitter_points(i+1,2), seg, iErr);
    emitter_seg_array = [emitter_seg_array; seg];
end
% add emmiter segments to object
%[Obj1, iErr] = Object_Create(Lorentz, Emitter, Obj1, iErr);
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
%[Obj2, iErr] = Object_Create(Lorentz, Collector, Obj2, iErr);
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
%[Obj3, iErr] = Object_Create(Lorentz, Gate, Obj3, iErr);
for i=1:length(gate_seg_array)
    [iErr] = Object_AddSegment(Lorentz, Gate, gate_seg_array(i), iErr);   
end %
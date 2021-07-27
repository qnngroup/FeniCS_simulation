clear;
clc;
close all;
%template with the right precision setup
template = 'C:\Users\John\Desktop\Lorentz_Scripted\template.tpl';

% Enhanced Fowler-Nordheim parameters
Work_Fun = 4.6;
beta = 100;
alpha = 1e-8;


% geometrical parameters
x_gap = 20;                 % gap between emitter and collector in the x direction 
y_gap = 100;                % gap between 2 gates in the y direction
w = 100;                    % width of emitter and collector
wg = 50;                    % width of gates
l = 200;                     % length of emitter
lg = 200;                   % length of gate
h = 50;                     % height of triangle
offset = y_gap/2+h+lg +3;      % offset in order to move things above the x axsis


%potentials
Ve = -5;                    %emitter potential
Vg = 10;                    %gate potential
Vc = 0;                     %collector potential

%variables initialization
iErr = int16(1);

Seg1 = int32(1);
Seg2 = int32(1);
Seg3 = int32(1);
Seg4 = int32(1);
Seg5 = int32(1);
Seg6 = int32(1);
Seg7 = int32(1);
Seg8 = int32(1);
Seg9 = int32(1);
Seg10 = int32(1);
Seg11 = int32(1);
Seg12 = int32(1);
Seg13 = int32(1);
Seg14 = int32(1);
Seg15 = int32(1);
Seg16 = int32(1);
Seg17 = int32(1);
Seg18 = int32(1);
Seg19 = int32(1);
Seg20 = int32(1);

Obj1 = int32(1);
Obj2 = int32(1);
Obj3 = int32(1);
Obj4 = int32(1);
Emitter = 'Emitter';
Emitter2 = 'Emitter2';
Gate_b = 'Gate_b';
Gate_t = 'Gate_t';
Collector = 'Collector';

numPointsReq = 100;
numPointsReq2 = 10;
numPointsAct = int32(1);
%numPointsAct2 = int32(1);


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
[iErr] = View_Set2DLimits(Lorentz, -5*w, -2*w, 5*w, 6*w, iErr);

%----------GEOMETRY-------------------------------------------------------------------------------------

% emitter geometry
[Seg1, iErr] = Geometry2D_CreateLine(Lorentz, -x_gap/2, offset, -x_gap/2-h, offset+w/2,  Seg1, iErr);
[Seg2, iErr] = Geometry2D_CreateLine(Lorentz, -x_gap/2-h, offset+ w/2,-x_gap/2-h-l, offset +w/2,  Seg2, iErr);
[Seg3, iErr] = Geometry2D_CreateLine(Lorentz, -x_gap/2-h-l, offset +w/2, -x_gap/2-h-l, offset -w/2,  Seg3, iErr);
[Seg4, iErr] = Geometry2D_CreateLine(Lorentz, -x_gap/2-h-l, offset -w/2, -x_gap/2-h, offset -w/2,  Seg4, iErr);
[Seg5, iErr] = Geometry2D_CreateLine(Lorentz, -x_gap/2-h, offset -w/2, -x_gap/2, offset, Seg5, iErr);
% object assignment
[Obj1, iErr] = Object_Create(Lorentz, Emitter, Obj1, iErr);
[iErr] = Object_AddSegment(Lorentz, Emitter, Seg1, iErr);
[iErr] = Object_AddSegment(Lorentz, Emitter, Seg2, iErr);
[iErr] = Object_AddSegment(Lorentz, Emitter, Seg3, iErr);
[iErr] = Object_AddSegment(Lorentz, Emitter, Seg4, iErr);
[iErr] = Object_AddSegment(Lorentz, Emitter, Seg5, iErr);

% Collector definition
[Seg6, iErr] = Geometry2D_CreateLine(Lorentz, x_gap/2, offset, x_gap/2+h, offset +w/2, Seg6, iErr);
[Seg7, iErr] = Geometry2D_CreateLine(Lorentz, x_gap/2+h, offset +w/2, x_gap/2+h+l, offset +w/2, Seg7, iErr);
[Seg8, iErr] = Geometry2D_CreateLine(Lorentz, x_gap/2+h+l, offset +w/2, x_gap/2+h+l, offset -w/2, Seg8, iErr);
[Seg9, iErr] = Geometry2D_CreateLine(Lorentz, x_gap/2+h+l, offset -w/2, x_gap/2+h, offset -w/2,  Seg9, iErr);
[Seg10, iErr] = Geometry2D_CreateLine(Lorentz, x_gap/2+h, offset -w/2, x_gap/2, offset, Seg10, iErr);

% Collector object definition
[Obj2, iErr] = Object_Create(Lorentz, Collector, Obj2, iErr);
[iErr] = Object_AddSegment(Lorentz, Collector, Seg6, iErr);
[iErr] = Object_AddSegment(Lorentz, Collector, Seg7, iErr);
[iErr] = Object_AddSegment(Lorentz, Collector, Seg8, iErr);
[iErr] = Object_AddSegment(Lorentz, Collector, Seg9, iErr);
[iErr] = Object_AddSegment(Lorentz, Collector, Seg10, iErr);

% bottom gate geometry
[Seg11, iErr] = Geometry2D_CreateLine(Lorentz, 0, offset -y_gap/2 , -wg/2, offset -y_gap/2-h, Seg11, iErr);
[Seg12, iErr] = Geometry2D_CreateLine(Lorentz, -wg/2, offset -y_gap/2-h, -wg/2, offset -y_gap/2-h-lg, Seg12, iErr);
[Seg13, iErr] = Geometry2D_CreateLine(Lorentz, -wg/2, offset -y_gap/2-h-lg, wg/2, offset -y_gap/2-h-lg, Seg13, iErr);
[Seg14, iErr] = Geometry2D_CreateLine(Lorentz, wg/2, offset -y_gap/2-h-lg, wg/2, offset -y_gap/2-h, Seg14, iErr);
[Seg15, iErr] = Geometry2D_CreateLine(Lorentz, wg/2, offset -y_gap/2-h, 0, offset -y_gap/2,  Seg15, iErr);
% Gate1 object definition
[Obj3, iErr] = Object_Create(Lorentz, Gate_b, Obj3, iErr);
[iErr] = Object_AddSegment(Lorentz, Gate_b, Seg11, iErr);
[iErr] = Object_AddSegment(Lorentz, Gate_b, Seg12, iErr);
[iErr] = Object_AddSegment(Lorentz, Gate_b, Seg13, iErr);
[iErr] = Object_AddSegment(Lorentz, Gate_b, Seg14, iErr);
[iErr] = Object_AddSegment(Lorentz, Gate_b, Seg15, iErr);

% top gate geometry
[Seg16, iErr] = Geometry2D_CreateLine(Lorentz, 0, offset +y_gap/2 , -wg/2, offset +y_gap/2+h, Seg16, iErr);
[Seg17, iErr] = Geometry2D_CreateLine(Lorentz, -wg/2, offset +y_gap/2+h, -wg/2, offset+ y_gap/2+h+lg, Seg17, iErr);
[Seg18, iErr] = Geometry2D_CreateLine(Lorentz, -wg/2, offset +y_gap/2+h+lg, wg/2, offset +y_gap/2+h+lg, Seg18, iErr);
[Seg19, iErr] = Geometry2D_CreateLine(Lorentz, wg/2, offset +y_gap/2+h+lg, wg/2, offset +y_gap/2+h, Seg19, iErr);
[Seg20, iErr] = Geometry2D_CreateLine(Lorentz, wg/2, offset +y_gap/2+h, 0, offset+ y_gap/2,  Seg20, iErr);
% Gate2 object definition
[Obj4, iErr] = Object_Create(Lorentz, Gate_b, Obj4, iErr);
[iErr] = Object_AddSegment(Lorentz, Gate_b, Seg16, iErr);
[iErr] = Object_AddSegment(Lorentz, Gate_b, Seg17, iErr);
[iErr] = Object_AddSegment(Lorentz, Gate_b, Seg18, iErr);
[iErr] = Object_AddSegment(Lorentz, Gate_b, Seg19, iErr);
[iErr] = Object_AddSegment(Lorentz, Gate_b, Seg20, iErr);

%----------ELECTRIC-------------------------------------------------------------------------------------

[iErr] = Model_SetAnalysisMode(Lorentz, 'Electric', iErr);


%set voltages
[iErr] = Physics_Set2DVoltage(Lorentz, Emitter, Ve, iErr);
[iErr] = Physics_Set2DVoltage(Lorentz, Gate_b, Vg, iErr);
[iErr] = Physics_Set2DVoltage(Lorentz, Gate_t, Vg, iErr);
[iErr] = Physics_Set2DVoltage(Lorentz, Collector, Vc, iErr);


%solve fields
[iErr] = Solution_RunSolver(Lorentz, iErr);

%----------TRAJECTORY-------------------------------------------------------------------------------------

[iErr] = Model_SetAnalysisMode(Lorentz, 'Trajectory', iErr);


%create emitter
[iErr] = Physics_CreateSurfaceEmitter(Lorentz, Emitter, iErr);
[iErr] = Physics_AddEmitterSurface_ByNumber(Lorentz, Emitter, Seg1, 1, iErr);
[iErr] = Physics_AddEmitterSurface_ByNumber(Lorentz, Emitter, Seg5, 1, iErr);
[numPointsAct,iErr] = Physics_SetNumberEmissionPoints(Lorentz, Emitter, numPointsReq, numPointsAct, iErr);


%create collectors
[iErr] = Physics_CreateSurfaceCollector(Lorentz, Collector, iErr);
[iErr] = Physics_AddCollectorSurface_ByNumber(Lorentz, Collector, Seg6, iErr);
[iErr] = Physics_AddCollectorSurface_ByNumber(Lorentz, Collector, Seg7, iErr);
[iErr] = Physics_AddCollectorSurface_ByNumber(Lorentz, Collector, Seg9, iErr);
[iErr] = Physics_AddCollectorSurface_ByNumber(Lorentz, Collector, Seg10, iErr);

[iErr] = Physics_CreateSurfaceCollector(Lorentz, Gate_b, iErr);
[iErr] = Physics_AddCollectorSurface_ByNumber(Lorentz, Gate_b, Seg11, iErr);
[iErr] = Physics_AddCollectorSurface_ByNumber(Lorentz, Gate_b, Seg12, iErr);
[iErr] = Physics_AddCollectorSurface_ByNumber(Lorentz, Gate_b, Seg14, iErr);
[iErr] = Physics_AddCollectorSurface_ByNumber(Lorentz, Gate_b, Seg15, iErr);

[iErr] = Physics_CreateSurfaceCollector(Lorentz, Gate_t, iErr);
[iErr] = Physics_AddCollectorSurface_ByNumber(Lorentz, Gate_t, Seg16, iErr);
[iErr] = Physics_AddCollectorSurface_ByNumber(Lorentz, Gate_t, Seg17, iErr);
[iErr] = Physics_AddCollectorSurface_ByNumber(Lorentz, Gate_t, Seg19, iErr);
[iErr] = Physics_AddCollectorSurface_ByNumber(Lorentz, Gate_t, Seg20, iErr);


%define emission properties
[numPointAct,iErr] = Physics_SetNumberEmissionPoints(Lorentz, Emitter, numPointsReq, numPointsAct, iErr);
[iErr] = Physics_SetEmissionRegimeTo_EnhancedFowlerNordheim(Lorentz, Emitter, Work_Fun, beta, alpha, iErr);


%solve trajectories
[iErr] = Solution_RunSolver(Lorentz, iErr);

Nrays=numPointsAct;
ray.X=zeros(Nrays,1);
ray.Y=zeros(Nrays,1);
ray.Z=zeros(Nrays,1);
ray.VX=zeros(Nrays,1);
ray.VY=zeros(Nrays,1);
ray.VZ=zeros(Nrays,1);
ray.E=zeros(Nrays,1);
ray.T=zeros(Nrays,1);
ray.D=zeros(Nrays,1);
ray.R=zeros(Nrays,1);
ray.C=zeros(Nrays,1);

X=zeros(Nrays,1);
Y=zeros(Nrays,1);
Z=zeros(Nrays,1);
VX=zeros(Nrays,1);
VY=zeros(Nrays,1);
VZ=zeros(Nrays,1);
E=zeros(Nrays,1);
T=zeros(Nrays,1);
D=zeros(Nrays,1);
R=zeros(Nrays,1);
C=zeros(Nrays,1);


%current evaluation

%first emission
BeamNumber =1;
EmissionNumber =1;
for i=1:numPointsAct
    flag=0;
    iPointNumber=1;
    while flag==0
    [VX(i), VY(i), VZ(i), E(i), iErr] = Analysis_GetKinematics_OnRay(Lorentz, BeamNumber, EmissionNumber, i, iPointNumber, VX(i), VY(i),VZ(i), E(i), iErr);
    [X(i), Y(i), Z(i), T(i), D(i), iErr] =   Analysis_GetCoordinateTimeDistance_OnRay(Lorentz, BeamNumber, EmissionNumber, i, iPointNumber, X(i), Y(i), Z(i), T(i), D(i), iErr);
    [R(i), C(i), iErr] = Analysis_GetRayRadiusAndCurrent(Lorentz, BeamNumber, EmissionNumber, i, R(i), C(i), iErr);
    if E(i)==0
        flag=1;
    else
        ray.VX(i)=VX(i); ray.VY(i)=VY(i); ray.VZ(i)=VZ(i); ray.E(i)=E(i); ray.X(i)=X(i); ray.Y(i)=Y(i); ray.Z(i)=Z(i); ray.T(i)=T(i); ray.D(i)=D(i); ray.R(i)=R(i); ray.C(i)=C(i);
    end
        iPointNumber=iPointNumber+1;
    end 
end


%X_Collector = x_eg+x_g+x_gc;
is_on6 = @(x, y) (abs(x*w/(2*h) - x_gap*w/(4*h) - y + offset) <= 0.05) & x<x_gap/2+h & x>x_gap/2;
is_on10 = @(x, y) (abs(-x*w/(2*h) + x_gap*w/(4*h) - y + offset) <= 0.05) & x<x_gap/2+h & x>x_gap/2;
is_on7 = @(x,y) (x>=x_gap/2+h) & (y==w/2+ offset) ;
is_on9 = @ (x,y) (x>=x_gap/2+h) & (y==-w/2 + offset);

Collector_vec = is_on6(ray.X, ray.Y) | is_on10(ray.X, ray.Y) | is_on7(ray.X, ray.Y) | is_on9(ray.X, ray.Y);
Collector_Current = -2*sum(ray.C(Collector_vec));
Gate_Current = -2*sum(ray.C(~Collector_vec));

Window_Close(Lorentz)

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
hg = 50;                    % height of triangle of gate
r_e = 5;                    % radius of emmiter tip nm
r_g = 10;                   % radius of gate tip
offsret = 0;      % offset in order to move things above the x axsis
r2d =180/pi;

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
Gate = 'Gate';
Collector = 'Collector';

numPointsReq = 100;
numPointsReq2 = 10;
numPointsAct = int32(1);
numPointsAct2 = int32(1);


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

%----------GEOMETRY-------------------------------------------------------------------------------------

% emitter geometry
theta = atan(0.5*w/h);
sin_theta = w/sqrt(w^2+4*h^2);
c_x = -x_gap/2 - r_e/sin_theta;
c_y = 0;
[Seg1, iErr] = Geometry2D_CreateLine(Lorentz, -x_gap/2-h-l, 0, c_x+r_e, 0, Seg1, iErr);
[Seg2, iErr] = Geometry2D_CreateArc_ByCenterRadiusAngles(Lorentz, c_x, c_y, r_e, 0, (pi/2-theta) * r2d, Seg2, iErr);
x1=0; y1=0;x2=0;y2=0;x3=0;y3=0;
[x1, y1, x2, y2, x3, y3, iErr] = Geometry2D_GetArcPointCoordinates(Lorentz, Seg2 , x1, y1, x2, y2, x3, y3, iErr);
[Seg3, iErr] = Geometry2D_CreateLine(Lorentz, x3, y3, -x_gap/2-h, w/2, Seg3, iErr);
[Seg4, iErr] = Geometry2D_CreateLine(Lorentz, -x_gap/2-h, w/2, -x_gap/2-h-l, w/2,  Seg4, iErr);
[Seg5, iErr] = Geometry2D_CreateLine(Lorentz, -x_gap/2-h-l, w/2, -x_gap/2-h-l, 0,  Seg5, iErr);



% object assignment
[Obj1, iErr] = Object_Create(Lorentz, Emitter, Obj1, iErr);
[iErr] = Object_AddSegment(Lorentz, Emitter, Seg2, iErr);
[iErr] = Object_AddSegment(Lorentz, Emitter, Seg3, iErr);
[iErr] = Object_AddSegment(Lorentz, Emitter, Seg4, iErr);
[iErr] = Object_AddSegment(Lorentz, Emitter, Seg5, iErr);

% Collector definition
c_x = x_gap/2 + r_e/sin_theta;
c_y = 0;
[Seg6, iErr] = Geometry2D_CreateLine(Lorentz, x_gap/2+l+h,  0, c_x-r_e, 0, Seg6, iErr);
[Seg7, iErr] = Geometry2D_CreateArc_ByCenterRadiusAngles(Lorentz, c_x, c_y, r_e, (pi/2+theta) * r2d, pi* r2d, Seg7, iErr);
x1=0; y1=0;x2=0;y2=0;x3=0;y3=0;
[x1c, y1c, x2, y2, x3, y3, iErr] = Geometry2D_GetArcPointCoordinates(Lorentz, Seg7 , x1, y1, x2, y2, x3, y3, iErr);
[Seg8, iErr] = Geometry2D_CreateLine(Lorentz, x1c, y1c, x_gap/2+h, w/2, Seg8, iErr);
[Seg9, iErr] = Geometry2D_CreateLine(Lorentz, x_gap/2+h, w/2, x_gap/2+h+l, w/2,  Seg9, iErr);
[Seg10, iErr] = Geometry2D_CreateLine(Lorentz, x_gap/2+h+l, w/2, x_gap/2+h+l, 0, Seg10, iErr);

% Collector object definition
[Obj2, iErr] = Object_Create(Lorentz, Collector, Obj2, iErr);
[iErr] = Object_AddSegment(Lorentz, Collector, Seg7, iErr);
[iErr] = Object_AddSegment(Lorentz, Collector, Seg8, iErr);
[iErr] = Object_AddSegment(Lorentz, Collector, Seg9, iErr);
[iErr] = Object_AddSegment(Lorentz, Collector, Seg10, iErr);

theta_g = atan(wg/(2*hg));
sin_theta_g = wg/sqrt(wg^2 +4*hg^2);
c_xg = 0;
c_yg = y_gap/2 + r_g/sin_theta_g;

% gate geometry
%Geometry2D_CreateLine(Lorentz, 0, y_gap/2+hg_lg , 0,c_yg -r_g,  Seg11, iErr);
[Seg11, iErr] = Geometry2D_CreateArc_ByCenterRadiusAngles(Lorentz, c_xg, c_yg, r_g, (pi+theta_g) * r2d, (3*pi/2+theta_g)* r2d, Seg11, iErr);
x1=0; y1=0;x2=0;y2=0;x3=0;y3=0;
[x1, y1, x2, y2, x3, y3, iErr] = Geometry2D_GetArcPointCoordinates(Lorentz, Seg11 , x1, y1, x2, y2, x3, y3, iErr);
[Seg12, iErr] = Geometry2D_CreateLine(Lorentz, x1, y1, -wg/2, y_gap/2+hg, Seg12, iErr);
[Seg13, iErr] = Geometry2D_CreateLine(Lorentz, -wg/2, y_gap/2+hg, -wg/2, y_gap/2+hg+lg, Seg13, iErr);
[Seg14, iErr] = Geometry2D_CreateLine(Lorentz, -wg/2, y_gap/2+hg+lg, wg/2, y_gap/2+hg+lg, Seg14, iErr);
[Seg15, iErr] = Geometry2D_CreateLine(Lorentz, wg/2, y_gap/2+hg+lg, wg/2, y_gap/2+hg,  Seg15, iErr);
[Seg16, iErr] = Geometry2D_CreateLine(Lorentz, wg/2, y_gap/2+hg, x3, y3, Seg16, iErr);
% Gate1 object definition
[Obj3, iErr] = Object_Create(Lorentz, Gate, Obj3, iErr);
[iErr] = Object_AddSegment(Lorentz, Gate, Seg11, iErr);
[iErr] = Object_AddSegment(Lorentz, Gate, Seg12, iErr);
[iErr] = Object_AddSegment(Lorentz, Gate, Seg13, iErr);
[iErr] = Object_AddSegment(Lorentz, Gate, Seg14, iErr);
[iErr] = Object_AddSegment(Lorentz, Gate, Seg15, iErr);
[iErr] = Object_AddSegment(Lorentz, Gate, Seg16, iErr);



%----------ELECTRIC-------------------------------------------------------------------------------------

[iErr] = Model_SetAnalysisMode(Lorentz, 'Electric', iErr);


%set voltages
[iErr] = Physics_Set2DVoltage(Lorentz, Emitter, Ve, iErr);
[iErr] = Physics_Set2DVoltage(Lorentz, Gate, Vg, iErr);
[iErr] = Physics_Set2DVoltage(Lorentz, Collector, Vc, iErr);


%solve fields
[iErr] = Solution_RunSolver(Lorentz, iErr);

%----------TRAJECTORY-------------------------------------------------------------------------------------

[iErr] = Model_SetAnalysisMode(Lorentz, 'Trajectory', iErr);


%create emitter
[iErr] = Physics_CreateSurfaceEmitter(Lorentz, Emitter, iErr);
[iErr] = Physics_AddEmitterSurface_ByNumber(Lorentz, Emitter, Seg2, 1, iErr);
[numPointsAct,iErr] = Physics_SetNumberEmissionPoints(Lorentz, Emitter, numPointsReq, numPointsAct, iErr);

[iErr] = Physics_CreateSurfaceEmitter(Lorentz, Emitter2, iErr);
[iErr] = Physics_AddEmitterSurface_ByNumber(Lorentz, Emitter2, Seg3, 1, iErr);
[numPointsAct2,iErr] = Physics_SetNumberEmissionPoints(Lorentz, Emitter2, numPointsReq2, numPointsAct2, iErr);


%create collectors
[iErr] = Physics_CreateSurfaceCollector(Lorentz, Collector, iErr);
[iErr] = Physics_AddCollectorSurface_ByNumber(Lorentz, Collector, Seg7, iErr);
[iErr] = Physics_AddCollectorSurface_ByNumber(Lorentz, Collector, Seg8, iErr);
[iErr] = Physics_AddCollectorSurface_ByNumber(Lorentz, Collector, Seg9, iErr);

[iErr] = Physics_CreateSurfaceCollector(Lorentz, Gate, iErr);
[iErr] = Physics_AddCollectorSurface_ByNumber(Lorentz, Gate, Seg11, iErr);
[iErr] = Physics_AddCollectorSurface_ByNumber(Lorentz, Gate, Seg12, iErr);
[iErr] = Physics_AddCollectorSurface_ByNumber(Lorentz, Gate, Seg13, iErr);
[iErr] = Physics_AddCollectorSurface_ByNumber(Lorentz, Gate, Seg15, iErr);
[iErr] = Physics_AddCollectorSurface_ByNumber(Lorentz, Gate, Seg16, iErr);


%define emission properties
%[numPointAct,iErr] = Physics_SetNumberEmissionPoints(Lorentz, Emitter, numPointsReq, numPointsAct, iErr);
[iErr] = Physics_SetEmissionRegimeTo_EnhancedFowlerNordheim(Lorentz, Emitter, Work_Fun, beta, alpha, iErr);
[iErr] = Physics_SetEmissionRegimeTo_EnhancedFowlerNordheim(Lorentz, Emitter2, Work_Fun, beta, alpha, iErr);


%solve trajectories
[iErr] = Solution_RunSolver(Lorentz, iErr);

Nrays=numPointsAct + numPointsAct2;
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

%second emission
BeamNumber =1;
EmissionNumber =2;
for j=numPointsAct+1:(numPointsAct+numPointsAct2)
    flag=0;
    iPointNumber=1;
    while flag==0
    [VX(j), VY(j), VZ(j), E(j), iErr] = Analysis_GetKinematics_OnRay(Lorentz, BeamNumber, EmissionNumber, j-numPointsAct, iPointNumber, VX(j), VY(j),VZ(j), E(j), iErr);
    [X(j), Y(j), Z(j), T(j), D(j), iErr] =   Analysis_GetCoordinateTimeDistance_OnRay(Lorentz, BeamNumber, EmissionNumber, j-numPointsAct, iPointNumber, X(j), Y(j), Z(j), T(j), D(j), iErr);
    [R(j), C(j), iErr] = Analysis_GetRayRadiusAndCurrent(Lorentz, BeamNumber, EmissionNumber, j-numPointsAct, R(j), C(j), iErr);
    if E(j)==0
        flag=1;
    else
        ray.VX(j)=VX(j); ray.VY(j)=VY(j); ray.VZ(j)=VZ(j); ray.E(j)=E(j); ray.X(j)=X(j); ray.Y(j)=Y(j); ray.Z(j)=Z(j); ray.T(j)=T(j); ray.D(j)=D(j); ray.R(j)=R(j); ray.C(j)=C(j);
    end
        iPointNumber=iPointNumber+1;
    end 
end

%X_Collector = x_eg+x_g+x_gc;
is_on6 = @(x, y) ((x - x1c).^2 + (y - y1c).^2 - r_g^2 <= 0.05) & x <= x1c; 
is_on7 = @(x, y) (abs(x*w/(2*h) - x_gap*w/(4*h) - y) <= 0.05) & x<x_gap/2+h & x>x1c;
is_on8 = @(x,y) (x>=x_gap/2+h) & (y==w/2 ) ;

Collector_vec = is_on6(ray.X, ray.Y) | is_on7(ray.X, ray.Y) | is_on8(ray.X, ray.Y);
Collector_Current = -2*sum(ray.C(Collector_vec));
Gate_Current = -2*sum(ray.C(~Collector_vec));

Solution_DeleteSolution(Lorentz, iErr);
Window_Close(Lorentz)

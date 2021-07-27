function [Collector_Current, Gate_Current] = pNVCT_Ic_Ig(Work_Fun, beta, alpha, x_e, x_t, x_eg, x_g, x_gc, x_c, r_e, r_g, y_g, w, theta, Ve, Vg, Vc, numPointsReq, numPointsReq2 )
    
%template with the right precision setup
template = 'C:\Users\John\Desktop\Lorentz_Scripted\template.tpl';

r2d =180/pi;


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

Obj1 = int32(1);
Obj2 = int32(1);
Obj3 = int32(1);
Emitter = 'Emitter';
Emitter2 = 'Emitter2';
Gate = 'Gate';
Collector = 'Collector';

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
[Seg1, iErr] = Geometry2D_CreateLine(Lorentz, -x_e, 0, -(r_e/sin(theta)-r_e), 0,  Seg1, iErr);
[Seg2, iErr] = Geometry2D_CreateArc_ByCenterRadiusAngles(Lorentz, -r_e/sin(theta), 0, r_e, 0, (pi/2-theta) * r2d, Seg2, iErr);
x1=0; y1=0;x2=0;y2=0;x3=0;y3=0;
[x1, y1, x2, y2, x3, y3, iErr] = Geometry2D_GetArcPointCoordinates(Lorentz, Seg2 , x1, y1, x2, y2, x3, y3, iErr);
[Seg3, iErr] = Geometry2D_CreateLine(Lorentz, x3, y3,-x_t, x_t*tan(theta),  Seg3, iErr);
[Seg4, iErr] = Geometry2D_CreateLine(Lorentz, -x_t, x_t*tan(theta),-x_e, x_e*tan(theta),  Seg4, iErr);
[Seg5, iErr] = Geometry2D_CreateLine(Lorentz,-x_e, x_e*tan(theta), -x_e, 0,  Seg5, iErr);
% emitter object definition
[Obj1, iErr] = Object_Create(Lorentz, Emitter, Obj1, iErr);
[iErr] = Object_AddSegment(Lorentz, Emitter, Seg2, iErr);
[iErr] = Object_AddSegment(Lorentz, Emitter, Seg3, iErr);
[iErr] = Object_AddSegment(Lorentz, Emitter, Seg4, iErr);
[iErr] = Object_AddSegment(Lorentz, Emitter, Seg5, iErr);


% gate geometry
[Seg6, iErr] = Geometry2D_CreateLine(Lorentz, x_eg, y_g+r_g , x_eg, w,  Seg6, iErr);
[Seg7, iErr] = Geometry2D_CreateLine(Lorentz, x_eg, w , x_eg + x_g, w,  Seg7, iErr);
[Seg8, iErr] = Geometry2D_CreateLine(Lorentz, x_eg+x_g, y_g+r_g , x_eg+x_g, w,  Seg8, iErr);
[Seg9, iErr] = Geometry2D_CreateLine(Lorentz, x_eg+r_g, y_g , x_eg+x_g-r_g, y_g,  Seg9, iErr);
[Seg10, iErr] = Geometry2D_CreateArc_ByCenterRadiusAngles(Lorentz, x_eg+r_g, y_g+r_g , r_g, pi * r2d, 3*pi/2 * r2d, Seg10, iErr);
[Seg11, iErr] = Geometry2D_CreateArc_ByCenterRadiusAngles(Lorentz, x_eg+x_g-r_g, y_g+r_g , r_g, 3*pi/2 * r2d, 0, Seg11, iErr);
% gate object definition
[Obj2, iErr] = Object_Create(Lorentz, Gate, Obj2, iErr);
[iErr] = Object_AddSegment(Lorentz, Gate, Seg6, iErr);
[iErr] = Object_AddSegment(Lorentz, Gate, Seg7, iErr);
[iErr] = Object_AddSegment(Lorentz, Gate, Seg8, iErr);
[iErr] = Object_AddSegment(Lorentz, Gate, Seg9, iErr);
[iErr] = Object_AddSegment(Lorentz, Gate, Seg10, iErr);
[iErr] = Object_AddSegment(Lorentz, Gate, Seg11, iErr);


% collector geometry
[Seg12, iErr] = Geometry2D_CreateLine(Lorentz, x_eg+x_g+x_gc, 0 , x_eg+x_g+x_gc, w,  Seg12, iErr);
[Seg13, iErr] = Geometry2D_CreateLine(Lorentz, x_eg+x_g+x_gc, w , x_eg+x_g+x_gc+x_c, w,  Seg13, iErr);
[Seg14, iErr] = Geometry2D_CreateLine(Lorentz, x_eg+x_g+x_gc+x_c, w , x_eg+x_g+x_gc+x_c, 0,  Seg14, iErr);
[Seg15, iErr] = Geometry2D_CreateLine(Lorentz, x_eg+x_g+x_gc, 0 , x_eg+x_g+x_gc+x_c, 0,  Seg15, iErr);
% collector object definition
[Obj3, iErr] = Object_Create(Lorentz, Collector, Obj3, iErr);
[iErr] = Object_AddSegment(Lorentz, Collector, Seg12, iErr);
[iErr] = Object_AddSegment(Lorentz, Collector, Seg13, iErr);
[iErr] = Object_AddSegment(Lorentz, Collector, Seg14, iErr);


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
[iErr] = Physics_AddCollectorSurface_ByNumber(Lorentz, Collector, Seg12, iErr);

[iErr] = Physics_CreateSurfaceCollector(Lorentz, Gate, iErr);
[iErr] = Physics_AddCollectorSurface_ByNumber(Lorentz, Gate, Seg6, iErr);
[iErr] = Physics_AddCollectorSurface_ByNumber(Lorentz, Gate, Seg7, iErr);
[iErr] = Physics_AddCollectorSurface_ByNumber(Lorentz, Gate, Seg8, iErr);
[iErr] = Physics_AddCollectorSurface_ByNumber(Lorentz, Gate, Seg9, iErr);
[iErr] = Physics_AddCollectorSurface_ByNumber(Lorentz, Gate, Seg10, iErr);
[iErr] = Physics_AddCollectorSurface_ByNumber(Lorentz, Gate, Seg11, iErr);


%define emission properties
[numPointAct,iErr] = Physics_SetNumberEmissionPoints(Lorentz, Emitter, numPointsReq, numPointsAct, iErr);
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

X_Collector = x_eg+x_g+x_gc;
Collector_Current = -2*sum(ray.C(ray.X==X_Collector));
Gate_Current = -2*sum(ray.C(ray.X~=X_Collector));

Window_Close(Lorentz)

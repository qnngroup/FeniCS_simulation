function [gate_current, collector_current] = getCurrents(Lorentz, Ve, Vc, Vg, Emitter, Emitter2, Gate, Collector, Work_Fun, beta, alpha, emitterSegs, gateSegs, collectorSegs, numPointsReq, numPointsReq2, is_collector, is_gate, iErr)

[iErr] = Model_SetAnalysisMode(Lorentz, 'Electric', iErr);


%set voltages
[iErr] = Physics_Set2DVoltage(Lorentz, Emitter, Ve, iErr);
[iErr] = Physics_Set2DVoltage(Lorentz, Gate, Vg, iErr);
[iErr] = Physics_Set2DVoltage(Lorentz, Collector, Vc, iErr);


%solve fields
[iErr] = Solution_RunSolver(Lorentz, iErr);

%----------TRAJECTORY-------------------------------------------------------------------------------------

[iErr] = Model_SetAnalysisMode(Lorentz, 'Trajectory', iErr);
numPointsAct = int32(1);
numPointsAct2 = int32(1);

Seg2 = emitterSegs(2);
%create emitter
[iErr] = Physics_CreateSurfaceEmitter(Lorentz, Emitter, iErr);
[iErr] = Physics_AddEmitterSurface_ByNumber(Lorentz, Emitter, Seg2, 1, iErr);
[numPointsAct,iErr] = Physics_SetNumberEmissionPoints(Lorentz, Emitter, numPointsReq, numPointsAct, iErr);

Seg3 = emitterSegs(3);
[iErr] = Physics_CreateSurfaceEmitter(Lorentz, Emitter2, iErr);
[iErr] = Physics_AddEmitterSurface_ByNumber(Lorentz, Emitter2, Seg3, 1, iErr);
[numPointsAct2,iErr] = Physics_SetNumberEmissionPoints(Lorentz, Emitter2, numPointsReq2, numPointsAct2, iErr);


%create collectors
Seg7 = collectorSegs(2);
Seg8 = collectorSegs(3);
Seg9 = collectorSegs(4);
[iErr] = Physics_CreateSurfaceCollector(Lorentz, Collector, iErr);
[iErr] = Physics_AddCollectorSurface_ByNumber(Lorentz, Collector, Seg7, iErr);
[iErr] = Physics_AddCollectorSurface_ByNumber(Lorentz, Collector, Seg8, iErr);
[iErr] = Physics_AddCollectorSurface_ByNumber(Lorentz, Collector, Seg9, iErr);

Seg11 = gateSegs(1);
Seg12 = gateSegs(2);
Seg13 = gateSegs(3);
Seg15 = gateSegs(5);
Seg16 = gateSegs(6);
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

Collector_vec = is_collector(ray.X, ray.Y);
collector_current = -2*sum(ray.C(Collector_vec));
gate_current = -2*sum(ray.C(is_gate(ray.X, ray.Y)));

Solution_DeleteSolution(Lorentz, iErr);
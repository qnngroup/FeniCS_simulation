function [gate_current, collector_current] = getCurrents_modified(Lorentz, Ve, Vc, Vg, Emitter, Gate, Collector, Work_Fun, beta, alpha, emitterSegs, gateSegs, collectorSegs, collector_points, gate_points, numPointsReq, iterNum, iErr)

[iErr] = Model_SetAnalysisMode(Lorentz, 'Electric', iErr);


%set voltages
[iErr] = Physics_Set2DVoltage(Lorentz, Emitter, Ve, iErr);
[iErr] = Physics_Set2DVoltage(Lorentz, Gate, Vg, iErr);
[iErr] = Physics_Set2DVoltage(Lorentz, Collector, Vc, iErr);


%solve fields
[iErr] = Solution_RunSolver(Lorentz, iErr);
% iterate over all xy in the problem
%[xReal(i), xImag(i), yReal(i), yImag(i), iErr] = Analysis_Get2DElectricField(Lorentz, X, Y, xReal(i), xImag(i), yReal(i), yImag(i), iErr)
%[xReals, xImags, yReals, yImags, iErr] = Analysis_Get2DElectricFields(Lorentz, Xarray, Yarray, xReals, xImags, yReals, yImags, iErr)

%
Xarr = linspace(-200, 200, 60);
Yarr = linspace(0, 300, 30);
%Zarr = zeros(1, 1);
%VoltR = zeros(100, 50);
%VoltI = zeros(100, 50);
%[VoltR, VoltI, iErr] = Analysis_GetVoltages(Lorentz, int32(length(Xarr)*length(Yarr)*length(Zarr)), Xarr, Yarr, Zarr, VoltR, VoltI, iErr);
%save(['Voltage_Ve_' num2str(Ve) 'Vc_' num2str(Vc) 'Vg_' num2str(Vg) 'it_0'], 'VoltR', 'VoltI') 

VoltR = zeros(length(Xarr), length(Yarr));
VoltI = zeros(length(Xarr), length(Yarr));
for xpt = 1:length(Xarr)
    for ypt = 1:length(Yarr)
        [VoltR(xpt, ypt), VoltI(xpt, ypt), iErr] = Analysis_Get2DVoltage(Lorentz, Xarr(xpt), Yarr(ypt), VoltR(xpt, ypt), VoltI(xpt, ypt), iErr);
    end 
end
save(['Data/Voltage_Ve_' num2str(Ve) 'Vc_' num2str(Vc) 'Vg_' num2str(Vg) 'iter_' num2str(iterNum) '.mat'], 'VoltR', 'VoltI') 
%----------TRAJECTORY-------------------------------------------------------------------------------------

[iErr] = Model_SetAnalysisMode(Lorentz, 'Trajectory', iErr);
numPointsActEm = ones(length(emitterSegs), 1, 'int32');
%numPointsActCol = ones(length(collectorSegs), 1, 'int32');
numPointsActGt = ones(length(gateSegs), 1, 'int32');
%numPointsAct2 = int32(1);

%emitter emission
for i=1:length(emitterSegs)
    segi = emitterSegs(i);
    [iErr] = Physics_CreateSurfaceEmitter(Lorentz, Emitter, iErr);
    [iErr] = Physics_AddEmitterSurface_ByNumber(Lorentz, Emitter, segi, 1, iErr);   
    [numPointsActEm(i),iErr] = Physics_SetNumberEmissionPoints(Lorentz, Emitter, numPointsReq, numPointsActEm(i), iErr);
end

% comment out collector emission.
% % Collector emission
% for i=1:length(collectorSegs)
%     segi = collectorSegs(i);
%     [iErr] = Physics_CreateSurfaceEmitter(Lorentz, Collector, iErr);
%     [iErr] = Physics_AddEmitterSurface_ByNumber(Lorentz, Collector, segi, 1, iErr);   
%     [numPointsActCol(i),iErr] = Physics_SetNumberEmissionPoints(Lorentz, Collector, numPointsReq, numPointsActCol(i), iErr);
% end

%create collectors
[iErr] = Physics_CreateSurfaceCollector(Lorentz, Collector, iErr);
for i = 1:length(collectorSegs)
    [iErr] = Physics_AddCollectorSurface_ByNumber(Lorentz, Collector, collectorSegs(i), iErr);
end 

% Gate emission
for i=1:length(gateSegs)
    segi = gateSegs(i);
    [iErr] = Physics_CreateSurfaceEmitter(Lorentz, Gate, iErr);
    [iErr] = Physics_AddEmitterSurface_ByNumber(Lorentz, Gate, segi, 1, iErr);   
    [numPointsActGt(i),iErr] = Physics_SetNumberEmissionPoints(Lorentz, Gate, numPointsReq, numPointsActGt(i), iErr);
end

[iErr] = Physics_CreateSurfaceCollector(Lorentz, Gate, iErr);
for i=1:length(gateSegs)
[iErr] = Physics_AddCollectorSurface_ByNumber(Lorentz, Gate, gateSegs(i), iErr);
end

%define emission properties
%[numPointAct,iErr] = Physics_SetNumberEmissionPoints(Lorentz, Emitter, numPointsReq, numPointsAct, iErr);
betaEm = beta;
alphaEm = alpha;
betaCol = beta/2;
alphaCol = beta/2;
betaGt = beta/2;
alphaGt = alpha/2;
[iErr] = Physics_SetEmissionRegimeTo_EnhancedFowlerNordheim(Lorentz, Emitter, Work_Fun, betaEm, alphaEm, iErr);
[iErr] = Physics_SetEmissionRegimeTo_EnhancedFowlerNordheim(Lorentz, Collector, Work_Fun, betaCol, alphaCol, iErr);
[iErr] = Physics_SetEmissionRegimeTo_EnhancedFowlerNordheim(Lorentz, Gate, Work_Fun, betaGt, alphaGt, iErr);


%solve trajectories
[iErr] = Solution_RunSolver(Lorentz, iErr);

%Nrays=numPointsAct + numPointsAct2;
Nrays = sum(numPointsActEm) + sum(numPointsActGt);
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

%emitter emission
BeamNumber =1;
EmissionNumber =1;
for i = 1:sum(numPointsActEm)
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
        
%Collector emission
BeamNumber =1;
EmissionNumber =2;
for j=sum(numPointsActEm)+1:sum(numPointsActEm)+ sum(numPointsActGt)
    flag=0;
    iPointNumber=1;
    while flag==0
    [VX(j), VY(j), VZ(j), E(j), iErr] = Analysis_GetKinematics_OnRay(Lorentz, BeamNumber, EmissionNumber, j-sum(numPointsActEm), iPointNumber, VX(j), VY(j),VZ(j), E(j), iErr);
    [X(j), Y(j), Z(j), T(j), D(j), iErr] =   Analysis_GetCoordinateTimeDistance_OnRay(Lorentz, BeamNumber, EmissionNumber, j-sum(numPointsActEm), iPointNumber, X(j), Y(j), Z(j), T(j), D(j), iErr);
    [R(j), C(j), iErr] = Analysis_GetRayRadiusAndCurrent(Lorentz, BeamNumber, EmissionNumber, j-sum(numPointsActEm), R(j), C(j), iErr);
    if E(j)==0
        flag=1;
    else
        ray.VX(j)=VX(j); ray.VY(j)=VY(j); ray.VZ(j)=VZ(j); ray.E(j)=E(j); ray.X(j)=X(j); ray.Y(j)=Y(j); ray.Z(j)=Z(j); ray.T(j)=T(j); ray.D(j)=D(j); ray.R(j)=R(j); ray.C(j)=C(j);
    end
        iPointNumber=iPointNumber+1;
    end 
end

% %Gate emission
% BeamNumber =1;
% EmissionNumber =2;
% for k=sum(numPointsActEm)+ sum(numPointsActCol)+1:sum(numPointsActEm)+ sum(numPointsActCol)+sum(numPointsActGt)
%     flag=0;
%     iPointNumber=1;
%     while flag==0
%     [VX(k), VY(k), VZ(k), E(k), iErr] = Analysis_GetKinematics_OnRay(Lorentz, BeamNumber, EmissionNumber, k-sum(numPointsActEm)-sum(numPointsActCol), iPointNumber, VX(k), VY(k),VZ(k), E(k), iErr);
%     [X(k), Y(k), Z(k), T(k), D(k), iErr] =   Analysis_GetCoordinateTimeDistance_OnRay(Lorentz, BeamNumber, EmissionNumber, k-sum(numPointsActEm)-sum(numPointsActCol), iPointNumber, X(k), Y(k), Z(k), T(k), D(k), iErr);
%     [R(k), C(k), iErr] = Analysis_GetRayRadiusAndCurrent(Lorentz, BeamNumber, EmissionNumber, k-sum(numPointsActEm)-sum(numPointsActCol), R(k), C(k), iErr);
%     if E(k)==0
%         flag=1;
%     else
%         ray.VX(k)=VX(k); ray.VY(k)=VY(k); ray.VZ(k)=VZ(k); ray.E(k)=E(k); ray.X(k)=X(k); ray.Y(k)=Y(k); ray.Z(k)=Z(k); ray.T(k)=T(k); ray.D(k)=D(k); ray.R(k)=R(k); ray.C(k)=C(k);
%     end
%         iPointNumber=iPointNumber+1;
%     end 
% end


collector_vec = zeros(size(ray.X), 'logical');

for i=1:length(collector_vec)
    collector_vec(i)= intersectsBoundary(collector_points, ray.X(i), ray.Y(i), 0.01);
end 

gate_vec = zeros(size(ray.X), 'logical');
for i=1:length(gate_vec)
    gate_vec(i)= intersectsBoundary(gate_points, ray.X(i), ray.Y(i), 0.01);
end 

collector_current = -2*sum(ray.C(collector_vec));
gate_current = -2*sum(ray.C(gate_vec));
X = ray.X;
Y = ray.Y;
C = ray.C;
save(['Data/Current_Ve_' num2str(Ve) 'Vc_' num2str(Vc) 'Vg_' num2str(Vg) 'iter_' num2str(iterNum) '.mat'], 'X', 'Y', 'C', 'collector_vec', 'gate_vec') 
Solution_DeleteSolution(Lorentz, iErr);
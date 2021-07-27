clear;
clc;
close all;

% Enhanced Fowler-Nordheim parameters
Work_Fun = 4.6;
beta = 100;
alpha = 1e-8;

% initial geometrical parameters
x_e = 100;                  %emitter size in x
x_t = 20;                   %emitter tip size in x
x_eg = 10;                  %emitter to gate gap size in x
x_g = 20;                   %gate size in x
x_gc = 20;                  %gate to collector gap size in x
x_c = 50;                   %collector size in x
r_e = 5;                    %emitter tip radius
r_g = 5;                    %gate corner radii
y_g = 10;                   %gate gap size in y
w = 100;                    %gate and collector size in y
theta = pi / 4;             %emitter tip half angle

%potentials
Ve = 0;                    %emitter potential
Vg = 10;                    %gate potential
Vc = 0;                     %collector potential\

numPointsReq = 500;
numPointsReq2 = 10;

Vg = linspace(0,10,10);
Vc = linspace(0,20,5);

for j=1:numel(Vg)
    for i=1:numel(Vc)
        [Collector_Current(i,j), Gate_Current(i,j)] = pNVCT_Ic_Ig(Work_Fun, beta, alpha, x_e, x_t, x_eg, x_g, x_gc, x_c, r_e, r_g, y_g, w, theta, Ve, Vg(j), Vc(i), numPointsReq, numPointsReq2 );
    end
end
plot(Vc,Collector_Current)

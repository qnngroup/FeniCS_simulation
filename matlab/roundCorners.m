function [new_coords, xs, ys] = roundCorners(coords, ind1, ind2, N, tol)

% doesn't account for lines that are purely horizontal or purely verticle;
x1 = coords(ind1, 1);
y1 = coords(ind1, 2);
x2 = coords(ind2, 1);
y2 = coords(ind2, 2);



if abs(y1 - coords(ind1+1, 2))>tol && abs(x1 - coords(ind1+1, 1)) > tol
    m1 = (y1 - coords(ind1+1, 2))/(x1 - coords(ind1+1, 1));
    m1_perp = -1/m1;
    b1 = y1 - m1_perp*x1;
    
end

if abs(y2 - coords(ind2-1, 2)) > tol && abs(x2 - coords(ind2-1, 1)) > tol 
    m2 = (y2 - coords(ind2-1, 2))/(x2 - coords(ind2-1, 1));
    m2_perp = -1/m2;
    b2 = y2 - m2_perp*x2;
end

if abs(y1 - coords(ind1+1, 2)) < tol && abs(x2 - coords(ind2-1, 1)) < tol
    xc = x1;
    yc = y2;
    if y2 < y1 
        phi_initial = pi/2;
        phi_final = pi;
    else
        phi_initial = -pi/2;
        phi_final = 0; 
    end
elseif abs(x1 - coords(ind1+1, 1)) < tol && abs(y2 - coords(ind2-1, 2)) < tol
    xc = x2;
    yc = y1;
    if y1 < y2
        phi_initial = 0;
        phi_final = pi/2;
    else 
        phi_initial = pi;
        phi_final = 3*pi/2;
    end  
elseif abs(y1 - coords(ind1+1, 2)) < tol
    xc = x1;
    yc = m2_perp*xc +b2;
    if y1 > yc
        phi_initial = pi/2;
        phi_final = atan((y2 - yc)/(x2 - xc)) +pi;
    else
        phi_initial = -pi/2;
        phi_final = atan((y2 - yc)/(x2 - xc));
    end
elseif abs(y2 -coords(ind2-1, 2)) < tol
    xc = x2;
    yc = m1_perp*xc +b1;
    if y2 > yc
        phi_initial = atan((y1 - yc)/(x1 - xc));
        phi_final = pi/2;
    else 
        phi_initial = atan((y1 - yc)/(x1 - xc)) +pi;
        phi_final = 3*pi/2;
    end
elseif abs(x1 - coords(ind1+1, 1)) < tol
    yc = y1;
    xc = (yc - b2)/m2_perp;
    if x1 > xc 
        phi_initial = 0;
    else
        phi_initial = pi;
    end
    if x2 > xc
        phi_final = atan((y2 - yc)/(x2 - xc));
    else
        phi_final = atan((y2 - yc)/(x2 - xc)) + pi;
    end
elseif abs(x2 - coords(ind2+1, 1)) < tol
    yc = y2;
    xc = (yc - b1)/m1_perp;
    if x1 > xc 
        phi_initial = atan((y1 - yc)/(x1 - xc)) ;
    else
        phi_initial = atan((y1 - yc)/(x1 - xc));
    end 
    if x2 > xc
        phi_final = 0;
    else
        phi_final = pi;
    end
else
    xc = (b1 - b2)/(m2_perp - m1_perp);
    yc = m1_perp*xc +b1;
    if (x1 - xc) < 0 
        phi_initial = atan((x1 - xc)/(y1 - yc)) + pi;
    else
        phi_initial = atan((y1 - yc)/(x1 - xc));
    end 

    if (x2 - xc) < 0 
        phi_final = atan((y2 - yc)/(x2 - xc)) + pi;
    else
        phi_final = atan((y2 - yc)/(x2 - xc));
    end 
end 

R = sqrt((xc-x1)^2 + (yc-y1)^2);
new_coords = zeros(N, 2);

if phi_final < phi_initial
    phi_final = phi_final +2*pi;
end

for i = 1:N
    phi_i = phi_initial + i*(phi_final - phi_initial)/(N+1);
    xi = xc + R*cos(phi_i);
    yi = yc + R*sin(phi_i);
    new_coords(i, :) = [xi, yi];
end 

phi = linspace(0, 2*pi);
xs = xc + R*cos(phi);
ys = yc + R*sin(phi);
    


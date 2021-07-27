function isOn = intersectsBoundary(points, x_star, y_star, tol)

isOn = false;
for i=1:length(points) -1
    if abs(points(i, 1) - points(i+1, 1)) < tol
        if abs(x_star - points(i, 1)) < tol && ((y_star <= points(i,2) && y_star >= points(i+1,2)) || (y_star >= points(i,2) && y_star <= points(i+1,2)))
            isOn = true;
        end 
    elseif abs(points(i, 2) - points(i+1, 2)) < tol
        if abs(y_star - points(i, 2)) < tol && ((x_star <= points(i,1) && x_star >= points(i+1,1)) || (x_star >= points(i,1) && x_star <= points(i+1,1)))
            isOn = true;
        end
    elseif (x_star <= points(i,1) && x_star >= points(i+1,1)) || (x_star >= points(i,1) && x_star <= points(i+1,1))
        m = (points(i,2) - points(i+1,2))/(points(i,1) - points(i+1,1));
        b = points(i,2) - m*points(i,1);
        if abs(y_star - m*x_star - b) < tol
            isOn = true;
        end 
    end
end
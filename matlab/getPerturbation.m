function pert = getPerturbation(array, index, R, tol)

if array(index, 2) == 0
   pert = [R*(rand - 0.5), 0];
% if not end or beginning point:
elseif index ~= 1 && index ~= length(array)
   if  abs(array(index+1,1) - array(index-1, 1)) < tol 
       pert = [0, R*(rand - 0.5)];
   elseif abs(array(index+1, 2) - array(index-1, 2)) < tol 
        pert = [R*(rand - 0.5), 0];
   else
       m = (array(index+1, 2) - array(index-1, 2))/(array(index+1,1) - array(index-1, 1));
       m_perp = -1/m;
       delx = R*(rand - 0.5);
       dely = delx*m_perp;
       pert = [delx, dely]/sqrt(1+m_perp^2);
   end
else 
    if  abs(array(2,1) - array(end-1, 1)) < tol
        pert = [0, R*(rand - 0.5)];
    elseif abs(array(2, 2) - array(end-1, 2)) < tol
        pert = [R*(rand - 0.5), 0];
    else
       m = (array(2, 2) - array(end-1, 2))/(array(2,1) - array(end-1, 1));
       m_perp = -1/m;
       delx = R*(rand - 0.5);
       dely = delx*m_perp;
       pert = [delx, dely]/sqrt(1+m_perp^2);
    end
end
        

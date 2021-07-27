function score = scoreGeom(Ve, Vg, Vc, collector_current, gate_current)

% Vce needs to be a column vactor 
Vce = (Vc - Ve)';
disp(size(Vce))
disp(size(Vg))
% Collector_current is size (len(Vc), len(Vg))
current_slope = collector_current(2:end, :) - collector_current(1:end-1, :);
for i=1:length(collector_current(1, :))
    current_slope(:, i) = Vce(1:end-1).*current_slope(:, i)./(Vce(2:end) - Vce(1:end-1));
end
disp(size(current_slope))
for j=1:length(current_slope(:, 1))
    current_slope(j, :) = current_slope(j, :).*Vg;
end

% Vg is a row vector
score = sum(current_slope(:)) + 0.5*sum(gate_current(:)) + 1/sum(collector_current(:)); 
function points = getSegmentationBetweenPoints(x1, y1, x2, y2, seglen)

Nsegs = round(sqrt((x1-x2)^2 +(y1-y2)^2)/seglen);
if Nsegs < 2
    points = [x1, y1; x2, y2];
else
    points = zeros(Nsegs, 2);
    points(1, :) = [x1, y1];
    %points(end, :) = [x2, y2];
    for i=1:Nsegs-1
        x_t = (1- i/Nsegs)*x1 + i/Nsegs*x2;
        y_t = (1- i/Nsegs) *y1 + i/Nsegs*y2;
        points(i+1, :) = [x_t, y_t];
    end 
end 
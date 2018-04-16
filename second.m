function out = second(result)
    dv = 0.5;
    minv = 5;
    maxv = 30;
    ddelta = 0.01;
    mindelta = -pi/3;
    maxdelta = pi/3;
    out = zeros(int32((maxv-minv)/dv), int32((maxdelta-mindelta)/ddelta));
    
    for i = 1:int32((maxv-minv)/dv)
        for j = 1:int32((maxdelta-mindelta)/ddelta)
            out(i, j) = result(i, j)>-1 && result(i, j)<1;
        end
    end
end


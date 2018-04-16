function [outbeta, outdelta, outerror] = sol()
    dbeta = 0.05;
    minbeta = -pi/3;
    maxbeta = pi/3;
    ddelta = 0.05;
    mindelta = -pi/3;
    maxdelta = pi/3;
    outbeta = zeros(int32((maxbeta-minbeta)/dbeta), 1);
    outdelta = zeros(int32((maxbeta-minbeta)/dbeta), 1);
    outerror = zeros(int32((maxbeta-minbeta)/dbeta), 1);
    options = optimoptions(@fmincon,'Algorithm', 'sqp-legacy', 'Display','none','OptimalityTolerance', 1e-5, 'StepTolerance', 1e-5, 'MaxFunctionEvaluations', 100, 'UseParallel', false);
 
    psidot = 17/50;
    v = 17;
    
   for i = 1:int32((maxbeta-minbeta)/dbeta)
       temp = zeros(1, int32((maxdelta-mindelta)/ddelta));
       beta = minbeta:dbeta:maxbeta;
    %    [pos, val] = fmincon(@(x)abs(func(psidot, beta, v(i), x)), 0, [], [], [], [], -pi/3, pi/3, [], options);
        
        delta = mindelta:ddelta:maxdelta;
        for j = 1:int32((maxdelta-mindelta)/ddelta)
            temp(j) = abs(func(psidot, beta(i), v, delta(j)));
        end
        [val, pos] = min(temp);
        outbeta(i) = beta(i);
        outdelta(i) = pos*ddelta;
        outerror(i) = val;
   end
end
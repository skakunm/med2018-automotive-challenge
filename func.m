function [out] = func(psidot, beta, v, delta)
    lf = 1.4;
    lr = 1.345;
    wf = 0.85;
    wr = 0.842;
    h = 0.35;
    m = 1300;
    Iz = 1400;
    
    r = v/psidot;    

    ax = v^2/r*sin(beta);
    ay = -v^2/r*cos(beta);
    
    dffy = ay*m*h*lr/((lf + lr)*2*wr);
    dfry = ay*m*h*lf/((lf + lr)*2*wf);
    dfx = ax*m*h/((lf + lr)*2);
    fflz = 3124 - dfx - dffy;
    ffrz = 3124 - dfx + dffy;
    frlz = 3252 + dfx - dfry;
    frrz = 3252 + dfx + dfry;
        
    realfl = (cross([0 0 psidot], [lf wf 0]) + [v*cos(beta) v*sin(beta) 0]) * angle2dcm(delta, 0, 0, 'ZYX');
    realfr = (cross([0 0 psidot], [lf -wf 0]) + [v*cos(beta) v*sin(beta) 0]) * angle2dcm(delta, 0, 0, 'ZYX');
    realrl = (cross([0 0 psidot], [-lr wr 0]) + [v*cos(beta) v*sin(beta) 0]);
    realrr = (cross([0 0 psidot], [-lr -wr 0]) + [v*cos(beta) v*sin(beta) 0]);
    
    afl = atan2(realfl(2), realfl(1));
    afr = atan2(realfr(2), realfr(1));
    arl = atan2(realrl(2), realrl(1));
    arr = atan2(realrr(2), realrr(1));
    
    options = optimoptions(@fmincon,'Algorithm', 'sqp-legacy', 'Display','none','OptimalityTolerance', 1e-12, 'StepTolerance', 1e-12, 'MaxFunctionEvaluations', 1000, 'UseParallel', false);
    kfl = fmincon(@(x)(abs(tiref(x, afl, fflz))), 0, [], [], [], [], -0.1, 0.1, [], options);
    kfr = fmincon(@(x)(abs(tiref(x, afl, fflz))), 0, [], [], [], [], -0.1, 0.1, [], options);
    
    [fflx, ffly] = tiref(kfl, afl, fflz);
    [ffrx, ffry] = tiref(kfr, -afr, ffrz);
    ffry = -ffry;
    
    ffx = fflx + ffrx;
    ffy = ffly + ffry;
    
    frx =  ffy*sin(delta) - ffx*cos(delta) - m*v*psidot*sin(beta);
    fry = -ffy*cos(delta) - ffx*sin(delta) + m*v*psidot*cos(beta);
    
%     [~, maxflx] = fminbnd(@(x)(-tirer(x, arl, frlz)), -1, 1);
%     [~, maxfrx] = fminbnd(@(x)(-tirer(x, -arr, frrz)), -1, 1);
%     [~, maxfly] = tirer(0, arl, frlz);
%     [~, maxfry] = tirer(0, -arr, frrz);
%     if abs(frx) > abs(maxflx) + abs(maxfrx)
%        out = abs(frx); 
%     elseif abs(fry) > abs(maxfly) + abs(maxfry)
%        out = abs(fry); 
%     else
        fopt = optimoptions('fsolve','Display','none');
        kr = fsolve(@(x)rearfunc(x, frx, fry, arl, arr, frlz, frrz), [0, 0], fopt);
        
        [frlx, frly] = tirer(kr(1), arl, frlz);
        [frrx, frry] = tirer(kr(2), -arr, frrz);
        frry = -frry;
        
        out = 1/Iz*(lf*((ffly + ffry)*cos(delta) + (fflx + ffrx)*sin(delta)) ...
        + wf*((ffly - ffry)*sin(delta) + (ffrx - fflx)*cos(delta)) ...
        + wr*(frrx - frlx) - lr*(frly + frry));
   % end
end


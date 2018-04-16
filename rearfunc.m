function f = rearfunc(x, frx, fry, arl, arr, frlz, frrz)
    krl = x(1);
    krr = x(2);
    [frlx, frly] = tirer(krl, arl, frlz);
    [frrx, frry] = tirer(krr, -arr, frrz);
    frry = -frry;
    
    f(1) = frx - frlx - frrx;
    f(2) = fry - frly - frry;
end


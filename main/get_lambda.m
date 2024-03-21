function lambda = get_lambda(Tf)
    f = 40000;
    g           = 1.4;
    R           = 8314/28.97;
    T           = (Tf-32)/1.8 + 273.15;
    a      = sqrt(1.4*R*T);
    lambda = a / f;
end
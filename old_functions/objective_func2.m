function point_pressure = objective_func2(phi, X, Y, Z, U, V, W, x, y, z)
    A=1; %signal amplitude
    P0=2.4; %transducer amplitude power
    a=0.01; %piston radius CHANGE LATER
    lambda=0.0086; %wavelength of sound in air (m)

    k=2*pi/lambda; %wavenumber
    
    b=[x-X; y-Y; z-Z];
    d=vecnorm(b,2,1);
    sine_func=vecnorm(cross(b, [U; V; W]),2,1)./d;
    Df=sinc(k*a*sine_func);
    P=P0*A*Df.*exp(1i*(phi+k*d))./d;
    point_pressure=abs(sum(P));
end
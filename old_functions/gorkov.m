function U = gorkov(complex_pressure)
    r=0.001; %radius of droplet (m)
    V=pi*r^2; %volume of particle
    c0=343; %speed of sound in air
    rho0=1.18; %density of air (Kg/m^3)
    cs=900; %speed of sound in droplet
    rhos=1000; %density of water
    omega=40000; %wave frequency

    K1 =(1/4)*V*((1/(c0^2*rho0))-(1/(cs^2*rhos)));
    K2 =(3/4)*V*((rho0-rhos)/(omega^2*rho0*(rho0+2*rhos)));
    
    P_mag=abs(complex_pressure);
    [px,py,pz]=gradient(complex_pressure);
    U = K1 * P_mag .^2 - K2 * (abs(px) .^2 + abs(py) .^2 + abs(pz) .^ 2);
end
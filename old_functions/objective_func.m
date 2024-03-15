function obj_func = objective_func(phases,transducers_x, transducers_y, transducers_z, U, V, W, points_x, points_y, points_z)
    w=0;
    s=10^10;
    complex_pressure=pressure_field(phases,transducers_x, transducers_y, transducers_z, U, V, W, points_x, points_y, points_z);
    gorkov_potentials=gorkov(complex_pressure);

    P_mag=abs(complex_pressure);
    obj_func=sum(w*P_mag-s*del2(gorkov_potentials),"all");
end
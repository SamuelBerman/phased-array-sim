function [transducers_x, transducers_y, transducers_z] = transducerSphereOld(bounds)
    [transducers_x, transducers_y, transducers_z]=sphere(7);
    transducers_x=transducers_x(:)*bounds{1}(2);
    transducers_y=transducers_y(:)*bounds{1}(2);
    transducers_z=transducers_z(:)*bounds{1}(2);
end
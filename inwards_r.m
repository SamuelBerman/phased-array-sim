function [U,V,W] = inwards_r(X,Y,Z)
    % must return unit vector
    vec_mag=vecnorm([X,Y,Z],2,2);
    U=-X./vec_mag;
    V=-Y./vec_mag;
    W=-Z./vec_mag;
end
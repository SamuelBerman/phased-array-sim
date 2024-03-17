function [U,V,W] = inwards_z(Z)
    % must return unit vector
    W=sign(Z)* -1;
    U=W*0;
    V=W*0;
end
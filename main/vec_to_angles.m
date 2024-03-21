function [theta_angle,phi_angle] = vec_to_angles(U,V,W)
    theta_angle=acos(W./sqrt(U.^2 +V.^2 + W.^2));
    phi_angle=atan2(V, U);
end
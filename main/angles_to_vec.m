function [U,V,W] = angles_to_vec(theta_angle,phi_angle)
    U=sin(theta_angle).*cos(phi_angle);
    V=sin(theta_angle).*sin(phi_angle);
    W=cos(theta_angle);
end
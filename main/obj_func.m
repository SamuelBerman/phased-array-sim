function loss = obj_func(phi, X, Y, Z, theta_angle, phi_angle, ox, oy, oz, lambda, a)
    phi = repmat(phi, 8, 1);
    theta_angle = repmat(theta_angle, 8, 1);
    phi_angle = repmat(phi_angle, 8, 1);
    X = [X; X; -X; -X; X; X; -X; -X];
    Y = [Y; -Y; Y; -Y; Y; -Y; Y; -Y];
    Z = [Z; Z; Z; Z; -Z; -Z; -Z; -Z];
    
    [U,V,W] = angles_to_vec(theta_angle,phi_angle);
    loss = sum(pressure_field_nogpu(phi, X, Y, Z, U, V, W, ox, oy, oz, lambda, a));
end
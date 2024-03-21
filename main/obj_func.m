function loss = obj_func(phases,tx, ty, tz, x_angle, y_angle, points_x, points_y, points_z, lambda)
    [U,V,W] = angles_to_vec(x_angle,y_angle);
    loss = -1*sum(pressure_field_nogpu(phases,tx, ty, tz, U, V, W, points_x, points_y, points_z, lambda));
end
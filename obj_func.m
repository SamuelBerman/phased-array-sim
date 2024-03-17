function loss = obj_func(phases,tx, ty, tz, U, V, W, points_x, points_y, points_z, lambda)
    loss = -1*sum(pressure_field(phases,tx, ty, tz, U, V, W, points_x, points_y, points_z, lambda));
end
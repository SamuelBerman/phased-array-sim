%clear all;

bounds={[-0.075 0.075],[-0.075 0.075],[-0.08 0.08]};
a = 0.005;

[X,Y,Z] = transducer_grid(16, 16, bounds, true);

[U,V,W] = inwards_z(Z);

[ox, oy, oz] = deal([0], [0], [0]);

[theta_angle, phi_angle] = vec_to_angles(U,V,W);

phi=X*0;

options = optimset('MaxFunEvals',100000,'MaxIter',10000,'Display','final','PlotFcns',@optimplotfval);
phi=fminunc(@(phases) obj_func(phases, X, Y, Z, theta_angle, phi_angle, ox, oy, oz, lambda, a), phi, options);
[U,V,W]=angles_to_vec(theta_angle,phi_angle);

%% COMPUTE PRESSURES

lambda = get_lambda(70);

slice_axes={[0 NaN NaN 1], [NaN NaN -0.03 1]};
slice_axes={};

render_bounds={[0 0.02],[-0.02 0.02],[-0.02 0.02]};
%render_bounds={[-0.1 0.1],[-0.1 0.1],[-0.1 0.1]};
%render_bounds={};

slice_data=render_slices(slice_axes, 500, bounds,phi,X,Y,Z,U,V,W, lambda, a);
volume_data=render_volumes(50, render_bounds, phi, X, Y, Z, U, V, W, lambda, a);

%% PLOTTING
plot_data(X,Y,Z,bounds,slice_axes,slice_data, render_bounds,volume_data,phi,U,V,W,ox,oy,oz, 200);

clear all;

lambda = get_lambda(70);
a=0.008;

bounds={[-0.1 0.1],[-0.1 0.1],[-0.1 0.1]};

[X,Y,Z] = plate_points(lambda*2, 2, [10 22]);

[U,V,W] = inwards_r(X,Y,Z);

[ox, oy, oz] = deal([0], [-0.02], [0]);

phi=X*0;

[theta_angle, phi_angle] = vec_to_angles(U,V,W);
options = optimset('MaxFunEvals',100000,'MaxIter',10000,'Display','final','PlotFcns',@optimplotfval);
%phi=fminunc(@(phases) obj_func(phases, X, Y, Z, theta_angle, phi_angle, ox, oy, oz, lambda, a), phi, options);

%% COMPUTE PRESSURES

slice_axes={[0 NaN NaN 1]};
render_bounds={[0 0.02],[-0.02 0.02],[-0.02 0.02]};
%render_bounds={[-0.1 0.1],[-0.1 0.1],[-0.1 0.1]};
render_bounds={};
%slice_axes={};
slice_bounds={[-0.1 0.1],[-0.1 0.1],[-0.08 0.08]};

slice_data=render_slices(slice_axes, 500, slice_bounds,phi,X,Y,Z,U,V,W, lambda, a);
volume_data=render_volumes(100, render_bounds, phi, X, Y, Z, U, V, W, lambda, a);

%% PLOTTING
plot_data(X,Y,Z,bounds,slice_axes,slice_data, render_bounds,volume_data,phi,U,V,W,ox,oy,oz);

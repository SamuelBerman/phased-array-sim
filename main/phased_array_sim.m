clear all;

lambda = get_lambda(70);

bounds={[-0.1 0.1],[-0.1 0.1],[-0.1 0.1]}; %set bounds of the work area

% define transducer locations
%[X,Y,Z] = deal([0 0],[0 0],[-0.1 0.1]); % define transducer locations manually
[X,Y,Z] = transducer_grid(2, 2, bounds); % grid layout
%[X,Y,Z] = transducer_sphere(3, 0.1); % sphere layout
%[X,Y,Z] = plate_points(lambda*2, 2, [10 22]); % place transducers in concentric rings

% transform locations
%[X,Y,Z] = translate(X, Y, Z, 0, 0, 0);
%[X,Y,Z] = rotate(X, Y, Z, 0, 0, 0);

% define transducer normal vectors
[U,V,W] = inwards_z(Z); % point up/down
%[U,V,W] = inwards_r(X,Y,Z); % point towards center

[ox, oy, oz] = deal([0], [0], [0]);

[theta_angle, phi_angle] = vec_to_angles(U,V,W);

phi=X*0;




%options = optimset('MaxFunEvals',100000,'MaxIter',10000,'Display','final','PlotFcns',@optimplotfval);
%options = optimset('MaxFunEvals',10000,'MaxIter',10000);
%phi=fminunc(@(phases) obj_func(phases, X, Y, Z, U, V, W, ox, oy, oz, lambda), phi, options);

% params_start = [X,Y,theta_angle,phi_angle];
% params_optimized=fminsearch(@(params) obj_func(phi, params(:,1), params(:,2), Z, params(:,3), params(:,4), ox, oy, oz, lambda), params_start, options);
% X=params_optimized(:,1);
% Y=params_optimized(:,2);
% theta_angle=params_optimized(:,3);
% phi_angle=params_optimized(:,4);


prob_func = @(X_,Y_,theta_angle_,phi_angle_) obj_func(phi, X_, Y_, Z, theta_angle_, phi_angle_, ox, oy, oz, lambda);
X_ = optimvar("X_");
Y_ = optimvar("Y_");
theta_angle_ = optimvar("theta_angle_");
phi_angle_ = optimvar("phi_angle_");
prob = optimproblem("Objective",prob_func(X_,Y_,theta_angle_,phi_angle_));
params_start.X_ = X;
params_start.Y_ = Y;
params_start.theta_angle_ = theta_angle;
params_start.phi_angle_ = phi_angle;
[solf,fvalf,eflagf,outputf] = solve(prob, params_start);



[U,V,W]=angles_to_vec(theta_angle,phi_angle);

%phi=mod(phi,2*pi);

%% COMPUTE PRESSURES

slice_axes={[0 NaN NaN 1]};
render_bounds={[0 0.02],[-0.02 0.02],[-0.02 0.02]};
%render_bounds={[-0.1 0.1],[-0.1 0.1],[-0.1 0.1]};
render_bounds={};
%slice_axes={};
slice_bounds={[-0.1 0.1],[-0.1 0.1],[-0.08 0.08]};

slice_data=render_slices(slice_axes, 1000, slice_bounds,phi,X,Y,Z,U,V,W, lambda);
volume_data=render_volumes(100, render_bounds, phi, X, Y, Z, U, V, W, lambda);

%% PLOTTING
plot_data(X,Y,Z,bounds,slice_axes,slice_data, render_bounds,volume_data,phi,U,V,W,ox,oy,oz);
%plot_image(slice_data{1}{4}(:,:,1));

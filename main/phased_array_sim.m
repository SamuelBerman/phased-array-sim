% compare to kenny's code
% better math/optimization from paper

clear all; % clear workspace of variables

f = 40000;
g           = 1.4;
R           = 8314/28.97;
Tf          = 70;
T           = (Tf-32)/1.8 + 273.15;
a      = sqrt(1.4*R*T);
lambda = a / f;


bounds={[-0.1 0.1],[-0.1 0.1],[-0.1 0.1]}; %set bounds of the work area

% define transducer locations
%[X,Y,Z] = deal([0 0],[0 0],[-0.1 0.1]); % define transducer locations manually
%[X,Y,Z] = transducer_grid(16, 16, bounds); % grid layout
[X,Y,Z] = transducer_sphere(3, 0.1); % sphere layout
%[X,Y,Z] = plate_points(lambda); % place transducers in concentric rings

% transform locations
%[X,Y,Z] = translate(X, Y, Z, 0, 0, 0);
%[X,Y,Z] = rotate(X, Y, Z, 0, 0, 0);

% define transducer normal vectors
%[U,V,W] = inwards_z(Z); % point up/down
[U,V,W] = inwards_r(X,Y,Z); % point towards center

[ox, oy, oz] = deal([0], [0], [0]);

phi=X*0;
options = optimset('MaxFunEvals',100000,'MaxIter',10000,'Display','final','PlotFcns',@optimplotfval);
%options = optimset('MaxFunEvals',10000,'MaxIter',10000);
%phi=fminunc(@(phases) obj_func(phases, X, Y, Z, U, V, W, ox, oy, oz, lambda), phi, options);
%phi=mod(phi,2*pi);

%% COMPUTE PRESSURES

slice_axes={[NaN NaN 0 1]};
render_bounds={[0 0.02],[-0.02 0.02],[-0.02 0.02]};
%slice_axes={};
render_bounds={};

slice_data=render_slices(slice_axes, 2000, bounds,phi,X,Y,Z,U,V,W, lambda);
volume_data=render_volumes(100, render_bounds, phi, X, Y, Z, U, V, W, lambda);

%% PLOTTING
%plot_data(X,Y,Z,bounds,slice_axes,slice_data, render_bounds,volume_data,phi,U,V,W,ox,oy,oz);
plot_image(slice_data{1}{4}(:,:,1));

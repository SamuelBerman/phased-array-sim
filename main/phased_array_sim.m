clear all;

lambda = get_lambda(70);
a=0.00001;

bounds={[-0.1 0.1],[-0.1 0.1],[-0.1 0.1]}; %set bounds of the work area

min_dist=0.025;
array_bounds={[-0.1+min_dist/2 0-min_dist/2],[-0.1+min_dist/2 0-min_dist/2],[-0.1 0]};

% define transducer locations
%[X,Y,Z] = deal([0 0],[0 0],[-0.1 0.1]); % define transducer locations manually
[X,Y,Z] = transducer_grid(3, 3, array_bounds, false); % grid layout
%[X,Y,Z] = transducer_sphere(3, 0.1); % sphere layout
%[X,Y,Z] = plate_points(lambda*2, 2, [10 22]); % place transducers in concentric rings

% transform locations
%[X,Y,Z] = translate(X, Y, Z, 0, 0, 0);
%[X,Y,Z] = rotate(X, Y, Z, 0, 0, 0);

% define transducer normal vectors
%[U,V,W] = inwards_z(Z); % point up/down
[U,V,W] = inwards_r(X,Y,Z); % point towards center

[ox, oy, oz] = deal([0 0], [-0.01 0.01], [0 0]);

[theta_angle, phi_angle] = vec_to_angles(U,V,W);

phi=X*0;

%options = optimset('MaxFunEvals',100000,'MaxIter',10000,'Display','final','PlotFcns',@optimplotfval);
%options = optimset('MaxFunEvals',10000,'MaxIter',10000);
%phi=fminunc(@(phases) obj_func(phases, X, Y, Z, U, V, W, ox, oy, oz, lambda), phi, options);


len=length(X);
prob_func = @(X_,Y_) obj_func(phi, X_, Y_, Z, theta_angle, phi_angle, ox, oy, oz, lambda, a);
X_ = optimvar("X_", len, 'LowerBound',min_dist/2, 'UpperBound',0.1-min_dist/2);
Y_ = optimvar("Y_", len, 'LowerBound',min_dist/2, 'UpperBound',0.1-min_dist/2);
prob = optimproblem("Objective", prob_func(X_,Y_), 'ObjectiveSense','maximize');
distances=optimconstr((len-1)*(len));

for t1=1:len
    for t2=1:len
        if t1~=t2
            distances(t1+(t2-1)*len) = sqrt((X_(t1)-X_(t2))^2 + (Y_(t1)-Y_(t2))^2) >= min_dist;
        end
    end
end

prob.Constraints.d_top = distances;
x0.X_=X;
x0.Y_=Y;
options = optimoptions('surrogateopt','Display','iter','PlotFcn','surrogateoptplot','UseParallel',true, 'MaxFunctionEvaluations',10000);
%options = optimoptions('fmincon','Display','iter','PlotFcn','optimplotx','MaxFunEvals',100000,'MaxIter',10000);
%options = optimoptions('patternsearch','Display','iter','PlotFcn',{@psplotbestf,@psplotfuncount});
%options = optimoptions('simulannealbnd','Display','iter','PlotFcn','saplotbestf');
%options = optimoptions("ga",'Display','iter',"PlotFcn","gaplotbestf");
%options = optimoptions("particleswarm",'Display','iter',"PlotFcn","pswplotbestf");

[solf,fvalf,eflagf,outputf] = solve(prob,x0,"Solver","surrogateopt",'Options',options,'ObjectiveDerivative',"finite-differences",'ConstraintDerivative',"finite-differences");

X=solf.X_;
Y=solf.Y_;

phi = repmat(phi, 8, 1);
theta_angle = [repmat(theta_angle, 4, 1); repmat(theta_angle, 4, 1)+pi];
phi_angle = repmat(phi_angle, 8, 1);
X = [X; X; -X; -X; X; X; -X; -X];
Y = [Y; -Y; Y; -Y; Y; -Y; Y; -Y];
Z = [Z; Z; Z; Z; -Z; -Z; -Z; -Z];


[U,V,W]=angles_to_vec(theta_angle,phi_angle);

%phi=mod(phi,2*pi);

%% COMPUTE PRESSURES

slice_axes={[0 NaN NaN 1]};
render_bounds={[0 0.02],[-0.02 0.02],[-0.02 0.02]};
%render_bounds={[-0.1 0.1],[-0.1 0.1],[-0.1 0.1]};
render_bounds={};
%slice_axes={};
slice_bounds={[-0.1 0.1],[-0.1 0.1],[-0.08 0.08]};

slice_data=render_slices(slice_axes, 200, slice_bounds,phi,X,Y,Z,U,V,W, lambda, a);
volume_data=render_volumes(100, render_bounds, phi, X, Y, Z, U, V, W, lambda, a);

%% PLOTTING
[U,V,W] = inwards_r(X,Y,Z);
plot_data(X,Y,Z,bounds,slice_axes,slice_data, render_bounds,volume_data,phi,U,V,W,ox,oy,oz);
%plot_image(slice_data{1}{4}(:,:,1));

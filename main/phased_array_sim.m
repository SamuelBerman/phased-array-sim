% TO-DO
% compare to kenny's code
% better math/optimization from paper
% graphing functions
% amplitude optimization?

clear all; % clear workspace of variables
%close all; % close all windows
%clc; % clear the console

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

slice_axes={[NaN NaN 0 0.5]};
render_bounds={[0 0.02],[-0.02 0.02],[-0.02 0.02]};
%slice_axes={};
%render_bounds={};

slice_data=render_slices(slice_axes, 400, bounds,phi,X,Y,Z,U,V,W, lambda);
volume_data=render_volumes(100, render_bounds, phi, X, Y, Z, U, V, W, lambda);

%% PLOTTING
f=figure('Name','Phased Array Sim','NumberTitle','off');
scatter3(X,Y,Z,200,phi*100,'filled','MarkerFaceAlpha',0.4);

xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
xlim(bounds{1})
ylim(bounds{2})
zlim(bounds{3})
hold on
%scatter3(ox,oy,oz,10,'black','filled','MarkerFaceAlpha',1);

light("Position",[0 0 1],"Style","infinite")

%scatter3(ox,oy,oz,100,'black');
quiver3(X,Y,Z,U,V,W,0.4,LineWidth=1,AutoScale="off",Color="black");
colormap(jet)

for i=1:numel(slice_axes) %log10 for log pressure
    s=slice(slice_data{i}{1},slice_data{i}{2},slice_data{i}{3},slice_data{i}{4},slice_axes{i}(1),slice_axes{i}(2),slice_axes{i}(3),'nearest');
    alpha(s,slice_axes{i}(4))
    material(s, "dull")
end

shading("interp")
cb=colorbar;
ylabel(cb,'Pressure [ Pa ] | Phase [ rad ]','FontSize',11,'Rotation',270)
hold off
daspect([1 1 1])
set(gcf,'Color',[1 1 1])

for i=1:numel(render_bounds)
    isoval=300;
    p1=patch(isosurface(volume_data{i}{1},volume_data{i}{2},volume_data{i}{3},volume_data{i}{4}, isoval), 'FaceAlpha', 1,'FaceColor', 'interp','LineStyle', 'none');
    patch(isocaps(volume_data{i}{1},volume_data{i}{2},volume_data{i}{3},volume_data{i}{4}, isoval), 'FaceAlpha', 1,'FaceColor', 'interp', 'LineStyle', 'none');
    isonormals(volume_data{i}{1},volume_data{i}{2},volume_data{i}{3},volume_data{i}{4},p1);
    isocolors(volume_data{i}{1},volume_data{i}{2},volume_data{i}{3},volume_data{i}{4},p1);

    light("Position",[-1 0 0],"Style","infinite")
    light("Position",[1 0 0],"Style","infinite")
    lighting gouraud;
end

% f=figure;
% colormap(hot);
% graph_slice=squeeze(U(:,2,:));
% I=imagesc(graph_slice,'Interpolation','bilinear');

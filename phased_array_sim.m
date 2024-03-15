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
global lambda;
lambda = a / f;


bounds={[-0.16 0.16],[-0.16 0.16],[-0.16 0.16]}; %set bounds of the work area

% define transducer locations
%[X,Y,Z] = deal([0 0],[0 0],[-0.1 0.1]); % define transducer locations manually
[X,Y,Z] = transducer_grid(16, 16, bounds); % grid layout
%[X,Y,Z] = transducer_sphere(3, 0.1); % sphere layout
%[X,Y,Z] = plate_points(); % place transducers in concentric rings

% transform locations
%[X,Y,Z] = translate(X, Y, Z, 0, 0, 0);
%[X,Y,Z] = rotate(X, Y, Z, 0, 0, 0);

% define transducer normal vectors
[U,V,W] = inwards_z(Z); % point up/down
%[U,V,W] = inwards_r(X,Y,Z); % point towards center

[ox, oy, oz] = deal([0 0 0 0 -0.01 0.01], [-0.01 0.01 0 0 0 0], [0 0 -0.01 0.01 0 0]);
%[ox, oy, oz] = transducer_sphere(5, 0.01);

phi=X*0;
options = optimset('MaxFunEvals',100000,'MaxIter',10000,'Display','final','PlotFcns',@optimplotfval);
%options = optimset('MaxFunEvals',10000,'MaxIter',10000);
phi=fminunc(@(phases) obj_func(phases, X, Y, Z, U, V, W, ox, oy, oz), phi, options);
%phi=mod(phi,2*pi);

%% COMPUTE PRESSURES

slice_axes={[0 NaN NaN 1], [NaN 0 NaN 1]};
render_bounds={[0 0.02],[-0.02 0.02],[-0.02 0.02]};
slice_axes={};
%render_bounds={};

slice_data=render_slices(slice_axes, 400, bounds,phi,X,Y,Z,U,V,W);
volume_data=render_volumes(100, render_bounds, phi, X, Y, Z, U, V, W);

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
    isoval=110;
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


function output = render_volumes(res, bounds, phi, X, Y, Z, U, V, W)
    output={};
    for i=1:numel(bounds)
        slice_range={linspace(bounds{1}(1),bounds{1}(2),res), linspace(bounds{2}(1),bounds{2}(2),res), linspace(bounds{3}(1),bounds{3}(2),res)};
        [x_grid,y_grid,z_grid]=meshgrid(slice_range{1},slice_range{2},slice_range{3});
        output{i}={x_grid, y_grid, z_grid, abs(gather(pressure_field(phi,X,Y,Z,U,V,W,x_grid,y_grid,z_grid)))};
    end
end


function output = render_slices(slice_axes, slice_res, bounds, phi, X, Y, Z, U, V, W)
    output={};
    for i=1:numel(slice_axes)
        for j=1:3
            slice_coord=slice_axes{i}(j);
            if ~isnan(slice_coord)
                slice_range={linspace(bounds{1}(1),bounds{1}(2),slice_res), linspace(bounds{2}(1),bounds{2}(2),slice_res), linspace(bounds{3}(1),bounds{3}(2),slice_res)};

                slice_range{j}=slice_coord;
                [x_grid,y_grid,z_grid]=meshgrid(slice_range{1},slice_range{2},slice_range{3});
                pressures=abs(pressure_field(phi,X,Y,Z,U,V,W,x_grid,y_grid,z_grid));

                slice_range{j}=linspace(slice_coord,slice_coord+0.0001,2);
                [x_grid,y_grid,z_grid]=meshgrid(slice_range{1},slice_range{2},slice_range{3});

                if j==1, dim=2; elseif j==2, dim=1; elseif j==3, dim=3; end
                output{i}={x_grid, y_grid, z_grid, cat(dim, pressures, pressures)};
            end
        end
    end
end

function loss = obj_func(phases,tx, ty, tz, U, V, W, points_x, points_y, points_z)
    loss = -1*sum(pressure_field(phases,tx, ty, tz, U, V, W, points_x, points_y, points_z));
end

function P_field = pressure_field(phases,tx, ty, tz, U, V, W, points_x, points_y, points_z)
    % define system parameters

    A=1; %signal amplitude
    P0=2.4; %transducer amplitude power
    a=0.016; %piston radius CHANGE LATER *pi?
    global lambda;
    k=2*pi/lambda; %wavenumber
    
    
    [p_i, t_i] = meshgrid(linspace(1, numel(points_x), numel(points_x)), linspace(1, numel(tx), numel(tx)));

    p_i=gpuArray(p_i);
    t_i=gpuArray(t_i);

    P_field=reshape(abs(sum(arrayfun(@point_pressure, p_i, t_i),1,"omitnan")), size(points_x, 1), size(points_x, 2), []);
    P_field=gather(P_field);

    function p = point_pressure(p_i, t_i)
        b_x=points_x(p_i)-tx(t_i);
        b_y=points_y(p_i)-ty(t_i);
        b_z=points_z(p_i)-tz(t_i);
        
        d=sqrt(b_x^2+b_y^2+b_z^2);

        cx=b_y*W(t_i)-b_z*V(t_i);
        cy=b_x*W(t_i)-b_z*U(t_i);
        cz=b_x*V(t_i)-b_y*U(t_i);

        sine_func=k*a*sqrt(cx^2+cy^2+cz^2)/d;
        
        Df=sin(sine_func)/(sine_func);
        
        p=P0*A*Df*exp(1i*(phases(t_i)+k*d))/d;
    end
end

function [tx, ty, tz] = transducer_grid(num_x, num_y, bounds)
    x=linspace(bounds{1}(1),bounds{1}(2),num_x);
    y=linspace(bounds{2}(1),bounds{2}(2),num_y);
    [tx,ty]=meshgrid(x,y);
    
    tz=bounds{3}(1)*ones(length(x),length(y));

    % tz=[tz -tz];
    % tx=[tx tx];
    % ty=[ty ty];
    
    tx=tx(:);
    ty=ty(:);
    tz=tz(:);
end



function [tx, ty, tz] = transducer_sphere(num_rings, r)
    [tx, ty, tz] = deal([],[],[]);
    
    start_theta=pi/12;
    d_theta=pi/12;
    phi_spacing=0.23*r;

    for ring_num=0:num_rings-1
        theta=start_theta+ring_num*d_theta;
        circ=2*pi*r*sin(theta);
        num_trans=floor(circ/phi_spacing);
        phis=linspace(0,2*pi,num_trans);
        
        x=r*sin(theta)*cos(phis);
        y=r*sin(theta)*sin(phis);
        z=repelem(r*cos(theta),num_trans);

        tx=[tx x];
        ty=[ty y];
        tz=[tz z];
    end

    tx=[tx tx];
    ty=[ty ty];
    tz=[tz -tz];
end

function [tx, ty, tz] = plate_points()
    [tx, ty, tz] = deal([],[],[]);

    num_rings=2;
    spacing=0.0165;
    z=0.1;
    global lambda;
    d_start=z+0.004;

    for ring_num=0:num_rings-1
        d=d_start+ring_num*lambda;
        r=sqrt(d^2-z^2)
        circ=2*pi*r;
        num_trans=floor(circ/spacing)
        thetas=linspace(0,2*pi,num_trans);
        
        x=r*cos(thetas);
        y=r*sin(thetas);

        tx=[tx x];
        ty=[ty y];
    end
    tz=repelem(z, numel(tx));

    tx=[tx tx];
    ty=[ty ty];
    tz=[tz -tz];
end

function [U,V,W] = inwards_z(Z)
    % must return unit vector
    W=sign(Z)* -1;
    U=W*0;
    V=W*0;
end

function [U,V,W] = inwards_r(X,Y,Z)
    % must return unit vector
    vec_mag=vecnorm([X,Y,Z],2,2);
    U=-X./vec_mag;
    V=-Y./vec_mag;
    W=-Z./vec_mag;
end

function [X, Y, Z] = rotate(Xi, Yi, Zi, tx, ty, tz)
    h=numel(Xi)/2;
    X=Xi;
    Y=Yi;
    Z=Zi;
    out = rotx(tx)*[Xi(1:h); Yi(1:h); Zi(1:h)];
    out = roty(ty)*out;
    out = rotz(tz)*out;
    X(1:h)=out(1,:);
    Y(1:h)=out(2,:);
    Z(1:h)=out(3,:);
end

function [X, Y, Z] = translate(Xi, Yi, Zi, tx, ty, tz)
    h=numel(Xi)/2;
    X=Xi;
    Y=Yi;
    Z=Zi;
    X(1:h)=Xi(1:h)+tx;
    Y(1:h)=Yi(1:h)+ty;
    Z(1:h)=Zi(1:h)+tz;
end

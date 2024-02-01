% TO-DO
% compare to kenny's code
% better math/optimization from paper
% graphing functions

clear all; %clear workspace of variables

bounds={[-0.1 0.1],[-0.1 0.1],[-0.1 0.1]}; %set bounds of the work area
render_bounds={[0 0.02],[-0.02 0.02],[-0.02 0.02]}; %set bounds of the render area

% define transducer locations
%[X,Y,Z] = transducer_grid(4, 4, bounds); % grid layout
%[X,Y,Z] = transducer_sphere_alt(); % sphere layout
%[X,Y,Z] = deal([0 0],[0 0],[-10*0.0086 10*0.0086]); % define transducer locations manually
%[X,Y,Z] = rotate(X,Y,Z, 0, 0, 0);
[X,Y,Z] = platePoints2();

% define transducer normal vectors
%[U,V,W] = inwards_z(Z); % point up/down
[U,V,W] = inwards_r(X,Y,Z); % point towards center
%[X,Y,Z] = translate(X, Y, Z, 0, 0, 0);

O = [0 0 0]; % point at which phases will be optimized

phi=X*0;
%phi(1)=pi;
options = optimset('MaxFunEvals',10000,'MaxIter',10000,'Display','final','PlotFcns',@optimplotfval);
%phi=fminsearch(@(phases) pressure_field(phases, X, Y, Z, U, V, W, O(1), O(2), O(3)), phi, options);


%% COMPUTE PRESSURES

slice_axes={};
%slice_axes={[0 NaN NaN 1]};
%output=render_slices(slice_axes, 1000, bounds,phi,X,Y,Z,U,V,W);

output={render_volume(100, render_bounds, phi, X, Y, Z, U, V, W)};


%% PLOTTING

f=figure('Name','Phased Array Sim','NumberTitle','off');
scatter3(X,Y,Z,200,phi,'filled','MarkerFaceAlpha',0.4);
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
xlim(bounds{1})
ylim(bounds{2})
zlim(bounds{3})
hold on
scatter3(O(1),O(2),O(3),100,'black');
quiver3(X,Y,Z,U,V,W,0.4,LineWidth=1,AutoScale="off",Color="black");
colormap(jet)

for i=1:numel(slice_axes)
    s=slice(output{i}{1},output{i}{2},output{i}{3},log10(output{i}{4}),slice_axes{i}(1),slice_axes{i}(2),slice_axes{i}(3),'nearest');
    alpha(s,slice_axes{i}(4))
end

shading("interp")
cb=colorbar;
ylabel(cb,'Log Pressure [ log10(Pa) ] | Phase [ rad ]','FontSize',11,'Rotation',270)
hold off
daspect([1 1 1])
set(gcf,'Color',[1 1 1])

isoval=300;
p1=patch(isosurface(output{1}{1},output{1}{2},output{1}{3},output{1}{4}, isoval), 'FaceAlpha', 1,'FaceColor', 'interp','LineStyle', 'none');
patch(isocaps(output{1}{1},output{1}{2},output{1}{3},output{1}{4}, isoval), 'FaceAlpha', 1,'FaceColor', 'interp', 'LineStyle', 'none');
camlight left;
camlight right;
lighting gouraud;
isonormals(output{1}{1},output{1}{2},output{1}{3},output{1}{4},p1);
isocolors(output{1}{1},output{1}{2},output{1}{3},output{1}{4},p1);

% f=figure;
% colormap(hot);
% graph_slice=squeeze(U(:,2,:));
% I=imagesc(graph_slice,'Interpolation','bilinear');


function output = render_volume(res, bounds, phi, X, Y, Z, U, V, W)
    slice_range={linspace(bounds{1}(1),bounds{1}(2),res), linspace(bounds{2}(1),bounds{2}(2),res), linspace(bounds{3}(1),bounds{3}(2),res)};
    [x_grid,y_grid,z_grid]=meshgrid(slice_range{1},slice_range{2},slice_range{3});
    output={x_grid, y_grid, z_grid, abs(pressure_field(phi,X,Y,Z,U,V,W,x_grid,y_grid,z_grid))};
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


function P_field = pressure_field(phases,transducers_x, transducers_y, transducers_z, U, V, W, points_x, points_y, points_z)
    % define system parameters

    A=1; %signal amplitude
    P0=2.4; %transducer amplitude power
    a=0.016*pi; %piston radius CHANGE LATER
    lambda=0.0086; %wavelength of sound in air (m)

    k=2*pi/lambda; %wavenumber
    
    [p_i, t_i] = meshgrid(linspace(1, numel(points_x), numel(points_x)), linspace(1, numel(transducers_x), numel(transducers_x)));

    p_i=gpuArray(p_i);
    t_i=gpuArray(t_i);

    P_field=reshape(abs(sum(arrayfun(@point_pressure, p_i, t_i),1)), size(points_x, 1), size(points_x, 2), []);

    function p = point_pressure(p_i, t_i)
        b_x=points_x(p_i)-transducers_x(t_i);
        b_y=points_y(p_i)-transducers_y(t_i);
        b_z=points_z(p_i)-transducers_z(t_i);
        
        d=sqrt(b_x^2+b_y^2+b_z^2);

        cx=b_y*W(t_i)-b_z*V(t_i);
        cy=b_x*W(t_i)-b_z*U(t_i);
        cz=b_x*V(t_i)-b_y*U(t_i);

        sine_func=k*a*sqrt(cx^2+cy^2+cz^2)/d;
        
        Df=sin(sine_func)/(sine_func);
        
        p=P0*A*Df*exp(1i*(phases(t_i)+k*d))/d;
    end
end


function point_pressure = objective_func2(phi, X, Y, Z, U, V, W, x, y, z)
    A=1; %signal amplitude
    P0=2.4; %transducer amplitude power
    a=0.01; %piston radius CHANGE LATER
    lambda=0.0086; %wavelength of sound in air (m)

    k=2*pi/lambda; %wavenumber
    
    b=[x-X; y-Y; z-Z];
    d=vecnorm(b,2,1);
    sine_func=vecnorm(cross(b, [U; V; W]),2,1)./d;
    Df=sinc(k*a*sine_func);
    P=P0*A*Df.*exp(1i*(phi+k*d))./d;
    point_pressure=abs(sum(P));
end


function U = gorkov(complex_pressure)
    r=0.001; %radius of droplet (m)
    V=pi*r^2; %volume of particle
    c0=343; %speed of sound in air
    rho0=1.18; %density of air (Kg/m^3)
    cs=900; %speed of sound in droplet
    rhos=1000; %density of water
    omega=40000; %wave frequency

    K1 =(1/4)*V*((1/(c0^2*rho0))-(1/(cs^2*rhos)));
    K2 =(3/4)*V*((rho0-rhos)/(omega^2*rho0*(rho0+2*rhos)));
    
    P_mag=abs(complex_pressure);
    [px,py,pz]=gradient(complex_pressure);
    U = K1 * P_mag .^2 - K2 * (abs(px) .^2 + abs(py) .^2 + abs(pz) .^ 2);
end


function obj_func = objective_func(phases,transducers_x, transducers_y, transducers_z, U, V, W, points_x, points_y, points_z)
    w=0;
    s=10^10;
    complex_pressure=pressure_field(phases,transducers_x, transducers_y, transducers_z, U, V, W, points_x, points_y, points_z);
    gorkov_potentials=gorkov(complex_pressure);

    P_mag=abs(complex_pressure);
    obj_func=sum(w*P_mag-s*del2(gorkov_potentials),"all");
end


function [transducers_x, transducers_y, transducers_z] = transducer_grid(num_x, num_y, bounds)
    x=linspace(bounds{1}(1),bounds{1}(2),num_x);
    y=linspace(bounds{2}(1),bounds{2}(2),num_y);
    [transducers_x,transducers_y]=meshgrid(x,y);
    transducers_x=horzcat(transducers_x,transducers_x);
    transducers_y=horzcat(transducers_y,transducers_y);
    z_bottom=bounds{3}(1)*ones(length(x),length(y));
    z_top=bounds{3}(2)*ones(length(x),length(y));
    transducers_z=horzcat(z_bottom,z_top);
    
    transducers_x=transducers_x(:);
    transducers_y=transducers_y(:);
    transducers_z=transducers_z(:);
end

function [transducers_x, transducers_y, transducers_z] = transducer_sphere(bounds)
    [transducers_x, transducers_y, transducers_z]=sphere(7);
    transducers_x=transducers_x(:)*bounds{1}(2);
    transducers_y=transducers_y(:)*bounds{1}(2);
    transducers_z=transducers_z(:)*bounds{1}(2);
end

function [tx, ty, tz] = transducer_sphere_alt()
    [tx, ty, tz] = deal([],[],[]);
    
    num_rings=3;
    r=0.1;
    start_theta=pi/12;
    d_theta=pi/12;
    phi_spacing=0.023;

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

function [tx, ty, tz] = platePoints()
    [tx, ty, tz] = deal([],[],[]);

    res=30;
    [x_vals, y_vals] = meshgrid(linspace(-0.1,0.1,res), linspace(-0.1,0.1,res));
    z=0.1;
    lambda=0.0086;

    for i=1:numel(x_vals)
        x=x_vals(i);
        y=y_vals(i);
        dist=sqrt(x^2+y^2+z^2);
        if mod(dist, lambda) < 0.001
            tx=[tx x];
            ty=[ty y];
        end
    end
    tz=repelem(z, numel(tx));

    tx=[tx tx];
    ty=[ty ty];
    tz=[tz -tz];
end

function [tx, ty, tz] = platePoints2()
    [tx, ty, tz] = deal([],[],[]);

    num_rings=3;
    spacing=0.02;
    z=0.1;
    lambda=0.0086;
    d_start=lambda*12.5;

    for ring_num=0:num_rings-1
        d=d_start+ring_num*lambda*1;
        r=sqrt(d^2-z^2);
        circ=2*pi*r;
        num_trans=floor(circ/spacing);
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

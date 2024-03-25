function P_field = pressure_field_nogpu(phases,tx, ty, tz, U, V, W, points_x, points_y, points_z, lambda, a)

    % define system parameters
    A=1; %signal amplitude
    P0=2.4; %transducer amplitude power
    k=2*pi/lambda; %wavenumber
    
    
    [p_i, t_i] = meshgrid(linspace(1, numel(points_x), numel(points_x)), linspace(1, numel(tx), numel(tx)));
    [p1,p2]=arrayfun(@point_pressure, p_i, t_i);
    p1=sum(p1,1);
    p2=sum(p2,1);
    abs_p=sqrt(p1.^2+p2.^2);
    P_field=reshape(abs_p, size(points_x, 1), size(points_x, 2), []);

    function [p1,p2] = point_pressure(p_i, t_i)
        b_x=points_x(p_i)-tx(t_i);
        b_y=points_y(p_i)-ty(t_i);
        b_z=points_z(p_i)-tz(t_i);
        
        d=sqrt(b_x^2+b_y^2+b_z^2);

        cx=b_y*W(t_i)-b_z*V(t_i);
        cy=b_x*W(t_i)-b_z*U(t_i);
        cz=b_x*V(t_i)-b_y*U(t_i);

        sine_func=k*a*sqrt(cx^2+cy^2+cz^2)/d;
        
        Df=sin(sine_func)/(sine_func);
        
        p1=P0*A*Df*cos(phases(t_i)+k*d) / d;
        p2=P0*A*Df*sin(phases(t_i)+k*d) / d;
    end
end
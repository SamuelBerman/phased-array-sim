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
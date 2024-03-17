function [tx, ty, tz] = plate_points(lambda)
    [tx, ty, tz] = deal([],[],[]);

    num_rings=2;
    spacing=0.0165;
    z=0.1;
    d_start=z+0.004;

    for ring_num=0:num_rings-1
        d=d_start+ring_num*lambda;
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
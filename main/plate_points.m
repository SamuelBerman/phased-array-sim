function [tx, ty, tz] = plate_points(lambda, num_rings, num_trans_arr)
    [tx, ty, tz] = deal([],[],[]);

    spacing=0.0165;
    z=0.1;
    d_start=z+0.004;

    for ring_num=0:num_rings-1
        d=d_start+ring_num*lambda;
        r=sqrt(d^2-z^2);
        circ=2*pi*r;
        if isempty(num_trans_arr)
            num_trans=floor(circ/spacing);
        else
            num_trans=num_trans_arr(ring_num+1);
        end
        thetas=linspace(0,2*pi,num_trans+1);
        
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
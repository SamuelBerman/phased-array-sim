function [tx, ty, tz] = platePointsOld()
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
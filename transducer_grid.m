function [tx, ty, tz] = transducer_grid(num_x, num_y, bounds)
    x=linspace(bounds{1}(1),bounds{1}(2),num_x);
    y=linspace(bounds{2}(1),bounds{2}(2),num_y);
    [tx,ty]=meshgrid(x,y);
    
    tz=bounds{3}(1)*ones(length(x),length(y));

    tz=[tz -tz];
    tx=[tx tx];
    ty=[ty ty];
    
    tx=tx(:);
    ty=ty(:);
    tz=tz(:);
end
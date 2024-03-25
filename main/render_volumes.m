function output = render_volumes(res, bounds, phi, X, Y, Z, U, V, W, lambda, a)
    output={};
    for i=1:numel(bounds)
        slice_range={linspace(bounds{1}(1),bounds{1}(2),res), linspace(bounds{2}(1),bounds{2}(2),res), linspace(bounds{3}(1),bounds{3}(2),res)};
        [x_grid,y_grid,z_grid]=meshgrid(slice_range{1},slice_range{2},slice_range{3});
        output{i}={x_grid, y_grid, z_grid, abs(gather(pressure_field(phi,X,Y,Z,U,V,W,x_grid,y_grid,z_grid, lambda, a)))};
    end
end
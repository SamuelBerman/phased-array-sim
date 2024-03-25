function output = render_slices(slice_axes, slice_res, bounds, phi, X, Y, Z, U, V, W, lambda, a)
    output={};
    for i=1:numel(slice_axes)
        for j=1:3
            slice_coord=slice_axes{i}(j);
            if ~isnan(slice_coord)
                slice_range={linspace(bounds{1}(1),bounds{1}(2),slice_res), linspace(bounds{2}(1),bounds{2}(2),slice_res), linspace(bounds{3}(1),bounds{3}(2),slice_res)};

                slice_range{j}=slice_coord;
                [x_grid,y_grid,z_grid]=meshgrid(slice_range{1},slice_range{2},slice_range{3});
                pressures=abs(pressure_field(phi,X,Y,Z,U,V,W,x_grid,y_grid,z_grid, lambda, a));

                slice_range{j}=linspace(slice_coord,slice_coord+0.0001,2);
                [x_grid,y_grid,z_grid]=meshgrid(slice_range{1},slice_range{2},slice_range{3});

                if j==1, dim=2; elseif j==2, dim=1; elseif j==3, dim=3; end
                output{i}={x_grid, y_grid, z_grid, cat(dim, pressures, pressures)};
            end
        end
    end
end
function plot_data(X,Y,Z,bounds,slice_axes,slice_data, render_bounds,volume_data,phi,U,V,W,ox,oy,oz, isoval)
    figure('Name','Phased Array Sim','NumberTitle','off');
    scatter3(X,Y,Z,200,phi*30,'filled','MarkerFaceAlpha',0.4);
    
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('z (m)')
    xlim(bounds{1})
    ylim(bounds{2})
    zlim(bounds{3})
    hold on
    
    scatter3(ox,oy,oz,100,'black');
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
        p1=patch(isosurface(volume_data{i}{1},volume_data{i}{2},volume_data{i}{3},volume_data{i}{4}, isoval), 'FaceAlpha', 1,'FaceColor', 'interp','LineStyle', 'none');
        patch(isocaps(volume_data{i}{1},volume_data{i}{2},volume_data{i}{3},volume_data{i}{4}, isoval), 'FaceAlpha', 1,'FaceColor', 'interp', 'LineStyle', 'none');
        isonormals(volume_data{i}{1},volume_data{i}{2},volume_data{i}{3},volume_data{i}{4},p1);
        isocolors(volume_data{i}{1},volume_data{i}{2},volume_data{i}{3},volume_data{i}{4},p1);
    
        light("Position",[-1 0 0],"Style","infinite")
        light("Position",[1 0 0],"Style","infinite")
        lighting gouraud;
    end
end
function plot_image(data)
    figure('Name','Phased Array Sim','NumberTitle','off');
    colormap(jet);
    imagesc(data,'Interpolation','bilinear');
    daspect([1 1 1])
    set(gcf,'Color',[1 1 1])
end
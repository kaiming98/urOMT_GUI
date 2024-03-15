function visualizeDisplacementVectors(tag)
%VISUALIZEPELINES visualize Lagrangian peclet-lines in matlab for data labelled by tag
configuration = Configuration3D(tag);
load(sprintf('%s/pathlines.mat',configuration.pathOutput));
if isempty(configuration.mask) 
    configuration.mask = Mask(configuration.pathMask,configuration.isMaskFilled,configuration.xRange,configuration.yRange,configuration.zRange);
    if configuration.do_resize
        configuration.mask = configuration.mask.resize(configuration.sizeFactor);
    end

    if configuration.dilate>0
        configuration.mask = configuration.mask.dilate(configuration.dilate,configuration.dilateDim);
    end
end
%% velocity flux vector 
figure,
strid = 10;
magnify = 1;
[x, y, z] = meshgrid(1:configuration.trueSize(2), 1:configuration.trueSize(1), 1:configuration.trueSize(3));
mskfv = isosurface(x,y,z,configuration.mask.contents,0.5);
mskp = patch(mskfv);
mskp.FaceColor = [.17,.17,.17];
mskp.FaceAlpha= 0.031;
mskp.EdgeColor = [.17,.17,.17];
mskp.EdgeAlpha= 0;
mskp.DisplayName = 'mask';

grid on, axis image
hold on,
q = quiver3(obj.startPointsInROI(1:strid:end,2),obj.startPointsInROI(1:strid:end,1),obj.startPointsInROI(1:strid:end,3),obj.displacementInROI(1:strid:end,2)*magnify,obj.displacementInROI(1:strid:end,1)*magnify,obj.displacementInROI(1:strid:end,3)*magnify,...
    'color','r','LineWidth',1.5,'MaxHeadSize',0.5,'AutoScale','off','DisplayName','flux vectors');
%title(sprintf('Velocity Flux Vectors'),'FontSize',20, 'Interpreter', 'none'), 
view([242.1011   14.4475])
ax = gca; ax.FontSize = 10; 
xlabel('x-axis','FontSize',18),ylabel('y-axis','FontSize',18),zlabel('z-axis','FontSize',18)
set(get(gca,'YLabel'),'Rotation',10);%xtickangle(20);
set(get(gca,'XLabel'),'Rotation',-30,'VerticalAlignment','middle');
set(gca,'Color',[1,1,1]), set(gcf,'unit','normalized','position',[0.4808 0.3259 0.2937 0.3625],'Color',[1,1,1]), set(gcf, 'InvertHardcopy', 'off')
grid on; axis image;
xticks(0:5:configuration.trueSize(1)); yticks(0:5:configuration.trueSize(2)); zticks(0:5:configuration.trueSize(3));
%// Compute the magnitude of the vectors
mags = sqrt(sum(cat(2, q.UData(:), q.VData(:), ...
            reshape(q.WData, numel(q.UData), [])).^2, 2));

%// Get the current colormap
currentColormap = colormap(jet);
[~, ~, ind] = histcounts(mags, size(currentColormap, 1));
cmap = uint8(ind2rgb(ind(:), currentColormap) * 255);
cmap(:,:,4) = 255;
cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);
set(q.Head, ...
    'ColorBinding', 'interpolated', ...
    'ColorData', reshape(cmap(1:3,:,:), [], 4).');   
set(q.Tail, ...
    'ColorBinding', 'interpolated', ...
    'ColorData', reshape(cmap(1:2,:,:), [], 4).');

view([242.1011   14.4475])
set(gca,'Color',[1,1,1]), set(gcf,'unit','normalized','position',[0.5152 0.2790 0.2943 0.3605],'Color',[1,1,1]), set(gcf, 'InvertHardcopy', 'off')
%colorbar('FontSize',18), 
grid on,
xlim([0 configuration.trueSize(1)]); ylim([0 configuration.trueSize(2)]); zlim([0 configuration.trueSize(3)])

cb = colorbar;
cb.Ticks = linspace(0, 800, 6);
cb.TickLabels = num2cell(0:160:800);
cb.FontSize = 8;
cb.Label.String = 'distance of movement (a.u.)';
title('velocity flux vector')

%saveas(gcf, sprintf('%s/%s/%s_LagFluxVector_E%02d_%02d.png',cfg.out_dir,cfg.outdir_Lag_v,cfg.tag,cfg.first_time,cfg.last_time+cfg.time_jump)); 
end


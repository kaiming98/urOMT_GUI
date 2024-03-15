function visualizePelines(tag)
%VISUALIZEPELINES visualize Lagrangian peclet-lines in matlab for data labelled by tag
configuration = Configuration3D(tag);
load(sprintf('%s/pathlines.mat',configuration.pathOutput));
 %% pelines
peSet = cell2mat(obj.PeInROI');
maxset = maxk(peSet(:),round(0.05*length(peSet(:))));
Pemax = maxset(end);
figure,
%Pemax = 1000;
Num = 100; % the larger, the more intervals in colors

SL2 = obj.positionInROI;
SL_Pe = obj.PeInROI;
nSL = length(SL2);
colors = jet(Num);
if isempty(configuration.mask) 
    configuration.mask = Mask(configuration.pathMask,configuration.isMaskFilled,configuration.xRange,configuration.yRange,configuration.zRange);
    if configuration.do_resize
        configuration.mask = configuration.mask.resize(configuration.sizeFactor);
    end

    if configuration.dilate>0
        configuration.mask = configuration.mask.dilate(configuration.dilate,configuration.dilateDim);
    end
end
for ind = 1:10:nSL
    SL_tmp = SL2{ind};
    SL_Pe_tmp = SL_Pe{ind};
    SL_Pe_tmp(SL_Pe_tmp>Pemax) = Pemax;
    SL_Pe_rk = round(SL_Pe_tmp/Pemax*Num);
    SL_Pe_rk(SL_Pe_rk<1) = 1;
    hlines = patch([SL_tmp(:,2);NaN],[SL_tmp(:,1);NaN],[SL_tmp(:,3);NaN],[1,1,1]);
    set(hlines,'FaceVertexCData',[colors(SL_Pe_rk,:);colors(SL_Pe_rk(end),:)],'EdgeColor','flat','FaceColor','none');
end
%
hold on;
[x, y, z] = meshgrid(1:configuration.trueSize(2), 1:configuration.trueSize(1), 1:configuration.trueSize(3));
mskfv = isosurface(x,y,z,configuration.mask.contents,0.5);
mskp = patch(mskfv);
mskp.FaceColor = [.37,.37,.37];
mskp.FaceAlpha= 0.1;
mskp.EdgeColor = [.37,.37,.37];
mskp.EdgeAlpha= 0;
view([242.1011   14.4475])
xticks(0:5:configuration.trueSize(1)); yticks(0:5:configuration.trueSize(2)); zticks(0:5:configuration.trueSize(3));
ax = gca; ax.FontSize = 10; 
xlabel('x-axis','FontSize',18),ylabel('y-axis','FontSize',18),zlabel('z-axis','FontSize',18)
set(get(gca,'YLabel'),'Rotation',10);%xtickangle(20);
set(get(gca,'XLabel'),'Rotation',-30,'VerticalAlignment','middle');
set(gca,'Color',[1,1,1]), set(gcf,'unit','normalized','position',[0.5152 0.2790 0.2943 0.3605],'Color',[1,1,1]), set(gcf, 'InvertHardcopy', 'off')
grid on; axis image;
xlim([0 configuration.trueSize(1)]); ylim([0 configuration.trueSize(2)]); zlim([0 configuration.trueSize(3)])
colormap('jet'); grid on;

cb = colorbar;
cb.Ticks = linspace(0, 1, 6);
cb.TickLabels = num2cell(round((0:Pemax/5:Pemax)));
cb.FontSize = 8;
cb.Label.String = '\it{Pe}';
%text(1,-3,25,'\it{Pe}','Rotation',90,'FontSize',25);
title('pelines')
xlim([0 configuration.trueSize(1)]); ylim([0 configuration.trueSize(2)]); zlim([0 configuration.trueSize(3)])
%saveas(gcf, sprintf('%s/%s/%s_LagPelines_E%02d_%02d.png',cfg.out_dir,cfg.outdir_Lag_v,cfg.tag,cfg.first_time,cfg.last_time+cfg.time_jump)); 
end


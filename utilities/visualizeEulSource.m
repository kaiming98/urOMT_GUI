function visualizeEulSource(tag)
%VISUALIZEEULSPEED show speed map in each time interval in matlab
configuration = Configuration3D(tag);
if ~configuration.addSource
    error('rOMT does not have relative source results');
end
load(sprintf('%s/eulerian.mat',configuration.pathOutput));
% Eul r
x = 1:configuration.trueSize(1);
y = 1:configuration.trueSize(2);
z = 1:configuration.trueSize(3);
src1 = obj.source(abs(obj.source)>0);
maxset = maxk(src1(:),round(0.05*length(src1(:))));
srcMax = maxset(end);
for i = 1:obj.nData
    figure,
    tmp = obj.source(:,:,:,i);
    tmp(tmp==0) = NaN;
    hs=slice(y,x,z,tmp,y,x,z);
    set(hs,'EdgeColor','none','FaceColor','interp','FaceAlpha',0.04);
    alpha('color'),alphamap([linspace(0,0.12,100)])
    grid off, box off, axis image
    xticks(0:10:configuration.trueSize(1)); yticks(0:5:configuration.trueSize(2)); zticks(0:5:configuration.trueSize(3));
    ax = gca; ax.FontSize = 10;
    xlabel('x-axis','FontSize',18),ylabel('y-axis','FontSize',18),zlabel('z-axis','FontSize',18)
    set(get(gca,'YLabel'),'Rotation',10);%xtickangle(20);
    set(get(gca,'XLabel'),'Rotation',-30,'VerticalAlignment','middle');
    clim([-srcMax,srcMax])
    colormap(bluewhitered_drdt)
    view([242.1011   14.4475])
    set(gca,'Color',[1,1,1]), set(gcf,'unit','normalized','position',[0.5152 0.2790 0.2943 0.3605],'Color',[1,1,1]), set(gcf, 'InvertHardcopy', 'off')
    grid on,
    xlim([0 configuration.trueSize(1)]); ylim([0 configuration.trueSize(2)]); zlim([0 configuration.trueSize(3)])

    cb = colorbar;
    cb.Ticks = linspace(-srcMax, srcMax, 5);
    cb.TickLabels = num2cell(round((-srcMax:srcMax/2:srcMax)*100)/100);
    cb.FontSize = 8;
    cb.Label.String = 'r (a.u.)';
    title(sprintf('Eulerian influx/clearance rate map E%d -> E%d',configuration.timeInitial+(i-1)*configuration.timeJump,configuration.timeInitial+i*configuration.timeJump));
    %text(1,-3,12,'relative source','Rotation',90,'FontSize',25);

    %saveas(gcf, sprintf('%s/%s/%s_EulAveR_E%02d_%02d.png',cfg.out_dir,cfg.outdir_Eul,cfg.tag,cfg.first_time,cfg.last_time+cfg.time_jump));
end
end
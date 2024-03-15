function visualizeEulSpeed(tag)
%VISUALIZEEULSPEED show speed map in each time interval in matlab
configuration = Configuration3D(tag);
load(sprintf('%s/eulerian.mat',configuration.pathOutput));
% Eul speed
x = 1:configuration.trueSize(1);
y = 1:configuration.trueSize(2);
z = 1:configuration.trueSize(3);
speed1 = obj.speed(obj.speed>0);
maxset = maxk(speed1(:),round(0.05*length(speed1(:))));
spdMax = maxset(end);
for i = 1:obj.nData
    figure,
    tmp = obj.speed(:,:,:,i);
    hs=slice(y,x,z,tmp,y,x,z);
    set(hs,'EdgeColor','none','FaceColor','interp','FaceAlpha',0.04);
    alpha('color'),alphamap(linspace(0,1,100))
    grid off, box off, axis image
    xticks(0:10:configuration.trueSize(1)); yticks(0:5:configuration.trueSize(2)); zticks(0:5:configuration.trueSize(3));
    ax = gca; ax.FontSize = 10;
    xlabel('x-axis','FontSize',18),ylabel('y-axis','FontSize',18),zlabel('z-axis','FontSize',18)
    set(get(gca,'YLabel'),'Rotation',10);%xtickangle(20);
    set(get(gca,'XLabel'),'Rotation',-30,'VerticalAlignment','middle');
    colormap(jet)
    clim([0,spdMax])
    view([242.1011   14.4475])
    set(gca,'Color',[1,1,1]), set(gcf,'unit','normalized','position',[0.5152 0.2790 0.2943 0.3605],'Color',[1,1,1]), set(gcf, 'InvertHardcopy', 'off')
    grid on,
    xlim([0 configuration.trueSize(1)]); ylim([0 configuration.trueSize(2)]); zlim([0 configuration.trueSize(3)])

    cb = colorbar;
    cb.Ticks = linspace(0, 1, 6);
    cb.TickLabels = num2cell(round((0:spdMax/5:spdMax)*100)/100);
    cb.FontSize = 8;
    cb.Label.String = 'speed (a.u.)';
    title(sprintf('Eulerian Speed map E%d -> E%d',configuration.timeInitial+(i-1)*configuration.timeJump,configuration.timeInitial+i*configuration.timeJump));
    %text(1,-3,15,'speed (a.u.)','Rotation',90,'FontSize',25);

    %saveas(gcf, sprintf('%s/%s/%s_EulAveSpeed_E%02d_%02d.png',cfg.out_dir,cfg.outdir_Eul,cfg.tag,cfg.first_time,cfg.last_time+cfg.time_jump));
end
end
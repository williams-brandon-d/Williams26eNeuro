function fig = plotPulseData(file,fieldName)

% norm_stim_ticks = [0 0.5 1];
plotLinewidth = 1;
% xlimits = [0 500]; % ms
xlimits = [file.pulse_start_time-100 file.pulse_stop_time+100]; % 100 ms before pulse to 100 ms after pulse

% switch file.info.comments
%     case 'DNQX after'
%         cweight = 1;
%     otherwise
%         cweight = 0.8;
% end

color = getCellcolor(file.info.cell_type,file.info.dataSet,file.info.comments);

% switch file.info.cell_type
%     case 'stellate'
%         color = cweight*[1 0 0];
%     case 'pyramidal'
%         color = cweight*[1 0.4 0];
%     case 'fast spiking'
%         color = cweight*[0 0 1];
% end

%     plotTitle = strjoin({info(ID_index).location info(ID_index).cell_type info(ID_index).cell_num };

raw_data = file.(fieldName).raw_data(:,file.(fieldName).trials_noArtifacts);
mean_raw_data = mean(raw_data,2);
mean_filt_data = mean(file.(fieldName).gamma_data(:,file.(fieldName).trials_noArtifacts),2);

file.(fieldName).trials_Artifacts = file.trials(file.(fieldName).spikes); % assumes 1 trial of data
raw_data_artifacts = file.(fieldName).raw_data(:,file.(fieldName).trials_Artifacts);

fontsize = loadFontSizes();

fig = figure('WindowState', 'maximized');
% fig = figure();

subplot(3,1,1)
% plot data without artifacts in gray
plot(file.time,raw_data,'Color',0.5*[1 1 1],'Linewidth',plotLinewidth);
% plot mean data in black
hold on;
han1 = plot(file.time,mean_raw_data,'Color','k','Linewidth',plotLinewidth,'DisplayName','Unfiltered');
% plot artifacts in pink
if ~isempty(raw_data_artifacts)
    plot(file.time,raw_data_artifacts,'Color','m','Linewidth',plotLinewidth);
end

if ~isempty(xlimits)
    xlim(xlimits); 
else
    xlim([min(file.time) max(file.time)]);
end

%    ylim([-455 110])

%     ylabel({'Membrane';'Current (pA)'},'FontSize',fontsize.ylabel,'FontWeight','bold')
%     title(sprintf('Trial-Averaged Data: nTrials = %d',nTrials),'FontSize',fontsize.tick)
%     title('Mean-subtracted Data','Fontsize',fontsize.tick)

ax2 = gca;
ax2.YAxis.FontSize = fontsize.tick;
ax2.XAxis.FontSize = fontsize.tick;
ax2.YAxis.FontWeight = 'bold';
ax2.XAxis.FontWeight = 'bold';
ax2.XAxis.Visible = 'off';
box off
set(gca,'visible','off')
% pbaspect([7 1 1])

% yscalebar subplot 1
yshift = 0.3;
xshift = 0.2;
Xlim = xlim;
Ylim = ylim;
yvalueshift = yshift*(Ylim(2)-Ylim(1));
xvalueshift = xshift*(Xlim(2)-Xlim(1));

y1scaleX = Xlim(1)*ones(1,2);
y1scaleY = Ylim(1)-yvalueshift + [0 round(0.4*range(mean_raw_data),1,'significant')];
hold on;
plot(y1scaleX,y1scaleY,'-k','Linewidth',2)

textX = y1scaleX(1) - 0.75*xvalueshift;
textY = y1scaleY(1) + range(y1scaleY)/2;
text(textX,textY,sprintf('%g %s',range(y1scaleY),file.(fieldName).data_units),'Fontsize',fontsize.ylabel,'Fontweight','bold')
hold off;

subplot(3,1,2)
displayName = sprintf('Filtered: %d-%d Hz',file.gamma_bandpass(1),file.gamma_bandpass(2));
han2 = plot(file.time, mean_filt_data,'-','Color',color,'Linewidth',plotLinewidth,'DisplayName',displayName);
if ~isempty(xlimits)
    xlim(xlimits); 
else
    xlim([min(file.time) max(file.time)]);
end

%    ylim([-132 60])
ax3 = gca;
ax3.YAxis.FontSize = fontsize.tick;
ax3.XAxis.FontSize = fontsize.tick;
ax3.YAxis.FontWeight = 'bold';
ax3.XAxis.FontWeight = 'bold';
%     ylabel({'Membrane';'Current (pA)'},'FontSize',fontsize.ylabel,'FontWeight','bold')
%     title(sprintf('Filtered: %d-%d Hz',data_bandpass(1),data_bandpass(2)),'FontSize',fontsize.tick)
box off
ax3.XAxis.Visible = 'off';
set(gca,'visible','off')

% time scalebar
yshift = 0.2;
Xlim = xlim;
Ylim = ylim;
yvalueshift = yshift*(Ylim(2)-Ylim(1));
tscaleX = [Xlim(1) Xlim(1)+100];
tscaleY = (Ylim(1)-yvalueshift)*ones(1,2);
hold on;
plot(tscaleX,tscaleY,'-k','Linewidth',2)

textX = tscaleX(1) + 0.1*range(tscaleX);
textY = tscaleY(1) - 0.9*yvalueshift;
text(textX,textY,sprintf('%g ms',range(tscaleX)),'Fontsize',fontsize.ylabel,'Fontweight','bold')
hold off;

%yscalebar subplot 2
xshift = 0.2;
Xlim = xlim;
Ylim = ylim;
xvalueshift = xshift*(Xlim(2)-Xlim(1));
y2scaleX = tscaleX(1)*ones(1,2);
%     y2scaleY = [tscaleY(1) tscaleY(1)+200];
y2scaleY = tscaleY(1) + [0 round(0.4*range(mean_filt_data),1,'significant')];
hold on;
plot(y2scaleX,y2scaleY,'-k','Linewidth',2)

textX = y2scaleX(1) - 0.75*xvalueshift;
textY = y2scaleY(1) + range(y2scaleY)/2;
text(textX,textY,sprintf('%g %s',range(y2scaleY),file.(fieldName).data_units),'Fontsize',fontsize.ylabel,'Fontweight','bold')

% fill space between time scale and y2 scale
%     fillscaleX = tscaleX(1);
%     plot(y2scaleX,y2scaleY,'-k','Linewidth',2)
hold off;

subplot(3,1,3)
han3 = area(file.time,file.stim.norm,'DisplayName','Light');
han3.FaceColor =  [91, 207, 244] / 255; % light blue
han3.EdgeColor = han3.FaceColor;
ylim([min(file.stim.norm) max(file.stim.norm)])

if ~isempty(xlimits)
    xlim(xlimits); 
else
    xlim([min(file.time) max(file.time)]);
end

ylim([0 1])
set(gca,'visible','off')
pbaspect([7 0.2 1])
pos = get(gca, 'Position'); % pos array = [x y width height]
pos(2) = 0.2;
set(gca, 'Position', pos)

%     xlabel('Time (ms)','FontSize',fontsize.xlabel,'FontWeight','bold')
%     ax1 = gca;
%     ax1.YTick = norm_stim_ticks;
%     ax1.YAxis.FontSize = fontsize.tick;
%     ax1.XAxis.FontSize = fontsize.tick;
%     ax1.YAxis.FontWeight = 'bold';
%     ax1.XAxis.FontWeight = 'bold';
%     ax1.XAxis.Visible = 'off';
%     ylabel({'Norm.';'Light Intensity'},'FontSize',fontsize.ylabel,'FontWeight','bold')
%     box off

%     ax1.Title.Visible = 'on';
%     title(plotTitle,'FontSize',fontsize.title,'FontWeight','bold');

legend(ax2,[han1 han2 han3],'Location','NorthOutside','Orientation','horizontal');
% lgd.FontSize = 14;
legend(ax2,'boxoff')

% pbaspect(ax2,[7 1 1])



end
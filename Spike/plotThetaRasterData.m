function figRaster = plotThetaRasterData(file)

% if strcmp(file.info.comments,{'DNQX after' ,'DNQX after 10 uM'})
% get ylim from DNQX before 

    if file.cycles(end) > 15
        plotCycles = 11:14;
    else
        plotCycles = 2:5;
    end
    
%     trial = 1;
    stim_ds = 10; % stim downsample for fill plotting

%     fontsize = 12;
%     plotLinewidth = 1;
% 
%     time_scaleBar = 100; % ms
%     tscaleshift = 25; % ms
% %     yscaleBar = 20; % mV
% 
%     switch dataType
%         case 'raw'
%             trial_data = file.(fieldName).raw_data(:,trial); 
%         case 'gamma'
%             trial_data = file.(fieldName).gamma_data(:,trial); 
%     end

trial = 1;

% plot params
LineFormat = struct();
LineFormat.LineWidth = 2;
LineFormat.LineStyle = '-';

% build spike times
time_ms = file.time*1000;
spike_times = {time_ms([file.spike_indices])'};

    stim_data = file.stim.raw_data(:,trial);
    stim_data = (stim_data - min(stim_data)) / range(stim_data); % 0 to 1
    stim_data = downsample(stim_data,stim_ds);

    time = file.time*1000; % ms
    stim_time = downsample(time,stim_ds);

    start_index = file.cycle_start_index_noArtifacts(plotCycles(1));
    stop_index = file.cycle_start_index_noArtifacts(plotCycles(end)+1);    
    xlimits = [time(start_index) time(stop_index)];
    

    [LineFormat.Color,~] = getCellcolor(file.info.cell_type,file.info.dataSet,file.info.comments);

    figRaster = figure;

    ax1 = subplot(2,1,1);
%     plot(time,trial_data,'Color',[color alpha],'Linewidth',plotLinewidth);
    [~, ~] = plotSpikeRaster(spike_times,'PlotType','vertline','LineFormat',LineFormat);

    if ~isempty(xlimits); xlim(xlimits); end % change x axis limits

    set(ax1,'visible','off')

%     % yscalebar subplot 1
%     Xlim = xlim;
% 
%     Ylim = ylim;
% %     Ylim = Ylim(1)+0.5*range(Ylim)+[-125 125]; % set DNQX before ylimits here
%     
%     y1scaleshift = 0.4*(Ylim(2)-Ylim(1));
% 
%     y1scaleX = Xlim(1)*ones(1,2);
%     yscaleBar = round(0.2*(Ylim(2)-Ylim(1)),1,'significant'); % 40 percent of range
%     y1scaleY = Ylim(1) + y1scaleshift + [0 yscaleBar];
% 
%     hold on;
%     plot(y1scaleX,y1scaleY,'-k','Linewidth',2);
%     drawnow;
% 
%     textX = y1scaleX(1);
%     textY = y1scaleY(1) + range(y1scaleY)/2;
%     y1scale = text(textX,textY,sprintf('%g %s ',range(y1scaleY),file.(fieldName).data_units));
%     y1scale.FontSize = fontsize;
%     y1scale.FontWeight = 'bold';
%     y1scale.HorizontalAlignment = 'right';
%     y1scale.VerticalAlignment = 'middle';
% 
%     % time scalebar
%     Xlim = xlim;
% %     Ylim = ylim;
% 
%     tscaleX = Xlim(1) + tscaleshift + [0 time_scaleBar];
%     tscaleY = Ylim(1)*ones(1,2) - 0.1*(Ylim(2)-Ylim(1));
% 
%     hold on;
%     plot(tscaleX,tscaleY,'-k','Linewidth',2); % plot time scalebar
%     drawnow;
% 
%     textX = tscaleX(1) + 0.5*time_scaleBar; % center of X scalebar
%     textY = tscaleY(1); % y position of X scalebar
%     xscale = text(textX,textY,sprintf('%g ms',range(tscaleX))); % plot text
%     xscale.FontSize = fontsize;
%     xscale.FontWeight = 'bold';
%     xscale.VerticalAlignment = 'top';
%     xscale.HorizontalAlignment = 'center';
%     hold off;
% 
%     ylim([tscaleY(1) Ylim(2)]);

    ax2 = subplot(2,1,2);
    han2 = fill(stim_time,stim_data,[91, 207, 244] / 255,'DisplayName','Light');
    han2.EdgeColor = 'none';
    han2.LineStyle = 'none';

    if ~isempty(xlimits); xlim(xlimits); end % change x axis limits
    ylim([0 1])

    set(ax2,'visible','off')
    drawnow;

    ax2pos = get(ax2, 'Position'); % pos array = [x y width height]
    ax2pos(2) = 0.45;
    ax2pos(4) = 0.05;
    set(ax2, 'Position', ax2pos)
    drawnow;


end
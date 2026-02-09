function fig = plotThetaVoltageGammaData(file,fieldName)

%     plotCycles = 9:12;
    plotCycles = 11:14;
    trial = 1;
    stim_ds = 10;

    fontsize = 12;
    plotLinewidth = 1;

    time_scaleBar = 100; % ms
    tscaleshift = 25; % ms

    yscaleBar = 20; % mV

    trial_data = file.(fieldName).gamma_data(:,trial); 
%     yscaleBar = round(yscalesize*range(trial_data),1,'significant'); % 40 percent of range

    stim_data = file.stim.raw_data(:,trial);
    stim_data = (stim_data - min(stim_data)) / range(stim_data); % 0 to 1
    stim_data = downsample(stim_data,stim_ds);

    time = file.time*1000; % ms
    stim_time = downsample(time,stim_ds);

    start_index = file.cycle_start_index_noArtifacts(plotCycles(1));
    stop_index = file.cycle_start_index_noArtifacts(plotCycles(end)+1);    
    xlimits = [time(start_index) time(stop_index)];
    
    color = getCellcolor(file.info.cell_type,file.info.dataSet,file.info.comments);

    fig = figure;

    ax1 = subplot(2,1,1);
    plot(time,trial_data,'Color',color,'Linewidth',plotLinewidth);

    if ~isempty(xlimits); xlim(xlimits); end % change x axis limits

    set(ax1,'visible','off')

    % yscalebar subplot 1
    Xlim = xlim;
    Ylim = ylim;
    y1scaleshift = 0.4*(Ylim(2)-Ylim(1));

    y1scaleX = Xlim(1)*ones(1,2);
%     y1scaleX = (Xlim(1)-0.1*(Xlim(2)-Xlim(1)))*ones(1,2);
    y1scaleY = Ylim(1) + y1scaleshift + [0 yscaleBar];

    hold on;
    plot(y1scaleX,y1scaleY,'-k','Linewidth',2)
    drawnow;
    
    textX = y1scaleX(1);
    textY = y1scaleY(1) + range(y1scaleY)/2;
    y1scale = text(textX,textY,sprintf('%g %s ',range(y1scaleY),file.(fieldName).data_units));
    y1scale.FontSize = fontsize;
    y1scale.FontWeight = 'bold';
    y1scale.HorizontalAlignment = 'right';
    y1scale.VerticalAlignment = 'middle';
    hold off;

    % time scalebar
    Xlim = xlim;
    Ylim = ylim;

    tscaleX = Xlim(1) + tscaleshift + [0 time_scaleBar];
    tscaleY = Ylim(1)*ones(1,2) - 0.1*(Ylim(2)-Ylim(1));

    hold on;
    plot(tscaleX,tscaleY,'-k','Linewidth',2); % plot time scalebar
    
    textX = tscaleX(1) + 0.5*time_scaleBar; % center of X scalebar
    textY = tscaleY(1); % y position of X scalebar
    xscale = text(textX,textY,sprintf('%g ms',range(tscaleX))); % plot text
    xscale.FontSize = fontsize;
    xscale.FontWeight = 'bold';
    xscale.VerticalAlignment = 'top';
    xscale.HorizontalAlignment = 'center';
    
    hold off;

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

%     ax1pos = get(ax1, 'Position'); % pos array = [x y width height]

%     fig.Units = "normalized";
%     figpos = get(fig,'Position');
%     figpos(2) = figpos(2) + ax2pos(2);
%     figpos(4) = ax1pos(2) + ax1pos(2);
%     set(fig,'Position',figpos);

end
function plotFIdata(data,plot_sweeps,dt,plotTitle,info)

    xlimits = [500 2500]; % Time (ms) for viewing only not analysis
    ylimits = [];
    plotTimeOffset = 1000;
    plotLinewidth = 1;
    tickfontsize = 15;
    labelfontsize = 20;
    titlefontsize = 25;

    nSamples = size(data,1);
    time = (0:nSamples-1)*dt*1000 - plotTimeOffset; % time in msec
    time = time';

    figure;
    plot(time,data(:,plot_sweeps),'-k','Linewidth',plotLinewidth)
    if ~isempty(xlimits); xlim(xlimits-plotTimeOffset); end
    if ~isempty(ylimits); ylim(ylimits); end
    xlabel('Time (ms)','FontSize',labelfontsize,'FontWeight','bold')
    ylabel('Membrane Voltage (mV)','FontSize',labelfontsize,'FontWeight','bold')
    ax = gca;
    ax.YAxis.FontSize = tickfontsize;
    ax.XAxis.FontSize = tickfontsize;
    ax.YAxis.FontWeight = 'bold';
    ax.XAxis.FontWeight = 'bold';
    box off
%         set(gca,'Visible','off') % turn everything off besides trace
    if isempty(plotTitle)
       han = title(strjoin({info.location info.cell_type info.cell_num }),'FontSize',titlefontsize,'FontWeight','bold');
    else
       han = title(plotTitle,'FontSize',titlefontsize,'FontWeight','bold');
    end
    set(han,'Visible','on') % make title visible

end
function fig = plotData(file,fieldName)

    cycle_start_index_Artifacts = file.cycle_start_index(file.(fieldName).spikes); % assumes 1 trial of data

    nSubplots = 4;
    hSubplots = gobjects(nSubplots,1);

    fig = figure('WindowState', 'maximized');

    for iSubplot =  1:nSubplots
    
        switch iSubplot 
            case 1
                data = file.(fieldName).raw_data(:,file.trials);
                color = 'k';
                Title = 'Raw';
            case 2
                data = file.(fieldName).gamma_data(:,file.trials);
                color = 'b';
                Title = 'Gamma';
            case 3 
                data = file.(fieldName).theta_data(:,file.trials);
                color = 'r';
                Title = 'Theta';
            case 4 
                data = file.stim.raw_data(:,file.trials);
                color = [91, 207, 244] / 255;
                Title = 'Stim';
        end

        if iSubplot == 4 
            yLabel = '';
        else
            yLabel = file.(fieldName).data_units;
        end

        hSubplots(iSubplot) = subplot(nSubplots,1,iSubplot);
        plot(file.time,data,'Color',color,'Linewidth',1);
        hold on
        title(Title,'FontSize',20,'FontWeight','bold')
        ylabel(yLabel);
%         xlim([-inf inf]);
        xlim([0 file.time(file.cycle_start_index(file.cycles(end))+file.cycle_length)]);
        ax = gca;
        ax.XAxis.FontSize = 15;
        ax.YAxis.FontSize = 15;
        ax.XAxis.FontWeight = 'bold';
        ax.YAxis.FontWeight = 'bold';
        box off
        
        for icycle = 1:numel(cycle_start_index_Artifacts)
            cycle_start = cycle_start_index_Artifacts(icycle);
            cycle_stop = cycle_start + file.cycle_length; 
            window_time = file.time(cycle_start:cycle_stop);
            window_data = data(cycle_start:cycle_stop);
            plot(window_time,window_data,'Color','m'); % highlight artifacts 
        end

    end

%     hSubplots(2) = subplot(4,1,2);
%     plot(file.time,file.(fieldName).gamma_data(:,file.trials),'b','Linewidth',1);
%     hold on
%     title('Gamma','FontSize',20,'FontWeight','bold');
%     ylabel(file.(fieldName).data_units);
%     xlim([-inf inf]);
%     ax = gca;
%     ax.XAxis.FontSize = 15;
%     ax.YAxis.FontSize = 15;
%     ax.XAxis.FontWeight = 'bold';
%     ax.YAxis.FontWeight = 'bold';
%     box off
% 
%     hSubplots(3) = subplot(4,1,3);
%     plot(file.time,file.(fieldName).theta_data(:,file.trials),'r','Linewidth',1);
%     hold on
%     title('Theta','FontSize',20,'FontWeight','bold');
%     ylabel(file.(fieldName).data_units);
%     xlim([-inf inf]);
%     ax = gca;
%     ax.XAxis.FontSize = 15;
%     ax.YAxis.FontSize = 15;
%     ax.XAxis.FontWeight = 'bold';
%     ax.YAxis.FontWeight = 'bold';
%     box off
% 
%     hSubplots(4) = subplot(4,1,4);
%     plot(file.time,file.stim.raw_data(:,file.trials),'Color',[91, 207, 244] / 255,'Linewidth',1);
%     hold on
%     title('Stim','FontSize',20,'FontWeight','bold');
%     xlim([-inf inf]);
%     ax = gca;
%     ax.XAxis.FontSize = 15;
%     ax.YAxis.FontSize = 15;
%     ax.XAxis.FontWeight = 'bold';
%     ax.YAxis.FontWeight = 'bold';
%     box off

    xlabel("Time (s)",'FontSize',20,'FontWeight','bold')
    % sgtitle(sprintf('CaMK2-ChR2 %g Hz Stim',lfp.stim_freq))

end
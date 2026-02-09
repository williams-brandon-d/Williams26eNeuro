function [file,figData,figRaster,figISIs] = detectPatchSpikes(file,fieldName,plotnum)

% not all spikes cross 0 mV threshold - changed to prominence threshold
% example PV Viral fast spiking 17 
% - lots of PV viral fast spiking cells have spikes < 0 mV 
% - could be due to pipette capacitance in bath setup on 2p scope
% fast spiking 45 

% detect spike params
% spike_threshold = 0; % mV
minPeakProm = 20; % mV
minPeakDist = 1; % ms
maxPeakWidth = 3; % ms

minPeakDistPoints = minPeakDist*(1e-3)*file.Fs;
maxPeakWidthPoints = maxPeakWidth*(1e-3)*file.Fs;

delta_phase = 2*pi/file.cycle_length;
cycle_phase = (-pi:delta_phase:pi)'; 

nCycles = numel(file.cycle_start_index_noArtifacts);

% plot params
LineFormat = struct();
LineFormat.LineWidth = 2;
LineFormat.LineStyle = '-';

% labelfontsize = 20;
tickfontsize = 15;
titlefontsize = 15;

[LineFormat.Color,~] = getCellcolor(file.info.cell_type,file.info.dataSet,file.info.comments); % RGBA

% find spike times 
for itrial = 1:file.nTrials
    trial = file.trials(itrial);
    trial_data = file.(fieldName).raw_data(:,trial);   

    [~,spike_indices] = findpeaks(trial_data,'MinPeakProminence',minPeakProm,'MinPeakDistance',minPeakDistPoints,'MaxPeakWidth',maxPeakWidthPoints);
    
    cycle_phase_cell = cell(nCycles,1);
    cycle_idx_cell = cell(nCycles,1);
    FRs = cell(nCycles,1);
    
    for icycle = 1:nCycles
%         cycle = file.cycles(icycle);
        cycle_indices = file.cycle_start_index_noArtifacts(icycle) + (0:file.cycle_length);
        cycle_data = trial_data(cycle_indices);
               
%         [~, cycle_spike_indices] = findpeaks(cycle_data,'MinPeakHeight',spike_threshold,'MinPeakDistance',minPeakDistPoints,'MaxPeakWidth',maxPeakWidthPoints);
        [~, cycle_spike_indices] = findpeaks(cycle_data,'MinPeakProminence',minPeakProm,'MinPeakDistance',minPeakDistPoints,'MaxPeakWidth',maxPeakWidthPoints);
    
        cycle_idx_cell{icycle,:} = cycle_spike_indices;
        cycle_phase_cell{icycle,:} = cycle_phase(cycle_spike_indices)';
        FRs{icycle,:} = 1./(diff(cycle_spike_indices)*file.dt); % Hz
    end
    
    emptyCells = cellfun('isempty',cycle_idx_cell);
    
    if plotnum
        % plot trial data with detected spike peaks
        plotTitle = strjoin({ file.info.cell_type file.info.cell_num sprintf('LED Input: %g mV',file.led_input) });

        figData = figure;
        plot(file.time,trial_data,'-','Color',LineFormat.Color)
        %             first_cycle_start = file.cycle_start_index_noArtifacts(1);
        %             xlim([file.time(first_cycle_start) file.time(first_cycle_start+3*file.cycle_length)])
        xlabel('Time (s)')
        ylabel('Membrane Voltage (mV)')
        title(plotTitle,'FontSize',titlefontsize,'FontWeight','bold')
        ax = gca;
        ax.YAxis.FontSize = tickfontsize;
        ax.YAxis.FontWeight = 'bold';
        ax.XAxis.FontSize = tickfontsize;
        ax.XAxis.FontWeight = 'bold';
        
        ylimits = ylim;
        ymax = ylimits(2);

        hold on; 
        for icycle = 1:nCycles
%             cycle = file.cycles(icycle);
            cycle_indices = file.cycle_start_index_noArtifacts(icycle) + (0:file.cycle_length);
            cycle_time = file.time(cycle_indices);
            spike_times = cycle_time(cycle_idx_cell{icycle});
            spike_y = repmat(ymax,numel(spike_times),1);
            plot(spike_times,spike_y,'|r')
        end
        hold off
        xlim([0 cycle_time(end)])




        figRaster = figure;
        if any(~emptyCells)
            [~, ~] = plotSpikeRaster(cycle_phase_cell,'PlotType','vertline','LineFormat',LineFormat);
        end
        xlim([cycle_phase(1) cycle_phase(end)])
        ylim([0 nCycles]+0.5)
        xticks([-pi -pi/2 0 pi/2 pi]);
        xticklabels({'-π','-π/2','0','π/2','π'});
        xlabel('Stim Theta Phase (rad)')
        ylabel('Theta Cycle #')
        title(plotTitle,'FontSize',titlefontsize,'FontWeight','bold')
        ax = gca;
        ax.YAxis.FontSize = tickfontsize;
        ax.YAxis.FontWeight = 'bold';
        ax.XAxis.FontSize = tickfontsize;
        ax.XAxis.FontWeight = 'bold';


%         edges = 0:2:50; % ms
        edges = 0:10:300; % Hz
        figISIs = figure;
%         histogram(cell2mat(FRs),edges,'FaceColor',LineFormat.Color);
        [counts,~] = histcounts(cell2mat(FRs),edges);
        delta_bin = edges(2) - edges(1);
        bin_midpoints = edges(1:end-1) + delta_bin/2;
        plot(bin_midpoints,counts/nCycles,'LineStyle','-','Color',LineFormat.Color,'Linewidth',1.5);
        xlabel('Interspike Firing Rate (Hz)')
        ylabel('Counts (1 / Theta Cycle)')
        title(plotTitle,'FontSize',titlefontsize,'FontWeight','bold')
        box off
        ax = gca;
        ax.YAxis.FontSize = tickfontsize;
        ax.YAxis.FontWeight = 'bold';
        ax.XAxis.FontSize = tickfontsize;
        ax.XAxis.FontWeight = 'bold';
    end
  
end

file.theta_spike_phases = cycle_phase_cell;
file.theta_spike_indices = cycle_idx_cell;
file.FRs = FRs;
file.spike_indices = spike_indices;


end
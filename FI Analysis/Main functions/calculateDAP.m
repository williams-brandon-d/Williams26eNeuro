function [average_dap,dap] = calculateDAP(data_sweep,spike_indices,Fs,pulse_start_index,pulse_stop_index,plotnum)

ISI_min = 20; % ms
stop_time = 6; % ms
stop_time_after_spike = stop_time/1000; % sec
xlimits = [500 2500];

nSpikes = length(spike_indices);
nSamples = length(data_sweep);
dt = 1/Fs;
time = (1:nSamples)*dt*1000;
%     dap = [];
dap = cell(nSpikes,1);
plot_flag = 0;

for iSpike = 1:nSpikes
    start_index = spike_indices(iSpike);
    
    if (start_index > pulse_stop_index) || (start_index < pulse_start_index)
        continue; % ignore spikes outside pulse
    end
    
    if iSpike < nSpikes
        ISI = (spike_indices(iSpike+1) - spike_indices(iSpike))*dt*1000;
        if ISI < ISI_min
            continue; % ignore small ISIs
        end
    end
    
    stop_index_after_spike = start_index + Fs*stop_time_after_spike;
    window_after_spike = data_sweep(start_index:stop_index_after_spike);
    
%         if iSpike == 1
%             dv = diff(window_after_spike);
%             dv2 = diff(dv);
%             figure;
%             plot(dv)
%             figure;
%             plot(dv2)
%         end
    
    [fAHP,fAHP_index] = min(window_after_spike);
    
    dap_window = window_after_spike(fAHP_index:end);
    [max_after_min,MaM_index] = max(dap_window);
    
    dap_iSpike = max_after_min - fAHP;
    
    if dap_iSpike > 30
       dap_iSpike = []; % ignore spikes ?
    end
    
%         dap = [dap dap_iSpike];
    dap{iSpike} = dap_iSpike;
    
    plot_flag = plot_flag + 1;
    
    if plotnum == 1 && plot_flag == 1
        figure;
        plot(time,data_sweep,'-k')
        hold on
        plot(time(start_index+fAHP_index-1),fAHP,'xb')
        plot( [time(1) time(end)],[fAHP fAHP],'--b' )
        plot(time(start_index+fAHP_index+MaM_index-2),max_after_min,'ob')
        plot( [time(1) time(end)],[max_after_min max_after_min],'--b' )
        plot( [time(start_index+fAHP_index+MaM_index-2) time(start_index+fAHP_index+MaM_index-2)],[fAHP max_after_min],'--k' )
        title( sprintf( 'dAP = %.1f mV', dap_iSpike ) )
        if ~isempty(xlimits); xlim(xlimits); end
        xlabel('Time (ms)')
        ylabel('Membrane Voltage (mV)')
    end
    
end
dap_mat = cell2mat(dap);
average_dap = mean(dap_mat);

end
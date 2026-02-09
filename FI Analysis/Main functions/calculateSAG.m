function [sag_hyper,sag_depol] = calculateSAG(data,pulse_start_index,pulse_stop_index,dt,plotnum)

    sag_sweep_index = 1; % 1st current sweep, -200 pA

%     pulse_start_sec = 1; % pulse starts at 1 sec
%     pulse_stop_sec = 2; % pulse ends at 2 sec    
%     pulse_start_index = Fs*pulse_start_sec + 1;
%     pulse_stop_index = Fs*pulse_stop_sec;

    data_hyper_sag = data( pulse_start_index:pulse_stop_index, sag_sweep_index );
    pulse_ss = median( data_hyper_sag ); % steady state of sweep during pulse
    [pulse_min,min_index] = min( data_hyper_sag ); % peak hyper sag
    sag_hyper = pulse_ss - pulse_min;
    
    depol_start_index = pulse_stop_index + 1;
    data_depol_sag = data( depol_start_index:end, sag_sweep_index );
    after_pulse_ss = median( data_depol_sag ); % steady state after pulse
    [after_pulse_max,~] = max( data_depol_sag ); % peak depol sag
    sag_depol = after_pulse_max - after_pulse_ss;
    
    if sag_depol > 50; sag_depol = []; end % ignore spike
    
    if plotnum == 1
        nSamples = size(data,1);
%         dt = 1/Fs;
        time = (1:nSamples)*dt*1000;
        time_hyper_peak = time( pulse_start_index + min_index );
%         time_depol_peak = time( depol_start_index + max_index );

        figure;
        plot(time, data(:,sag_sweep_index),'-k','Linewidth', 2)
        hold on
        plot([time(1) time(end)],[pulse_ss pulse_ss],'--r')
        plot([time(1) time(end)],[pulse_min pulse_min],'--r')
        plot([time_hyper_peak time_hyper_peak],[pulse_min pulse_ss],'--r')
        plot(time_hyper_peak,pulse_min,'or')
%         if ~isempty( sag_depol )
%             plot([time(1) time(end)],[after_pulse_ss after_pulse_ss],'--b')
%             plot([time(1) time(end)],[after_pulse_max after_pulse_max],'--b')
%             plot([time_depol_peak time_depol_peak],[after_pulse_max after_pulse_ss],'--b')
%             plot(time_depol_peak,after_pulse_max,'ob')
%         end
        xlabel('Time (ms)')
        ylabel('Membrane Voltage (mV)')
        title( sprintf('Hyper. SAG = %.1f mV',sag_hyper) )
%         title( sprintf('Hyper. SAG = %.1f mV, Depol. SAG = %.1f mV',sag_hyper,sag_depol) )
    end
    
end
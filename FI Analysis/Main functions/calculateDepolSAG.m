function sag_depol = calculateDepolSAG(data_sweep,Fs,pulse_start_index,pulse_stop_index,plotnum)

    depol_stop_sec = 1.2;
    depol_stop_index = uint64(Fs*depol_stop_sec);
    
    data_pulse = data_sweep( pulse_start_index:pulse_stop_index );
    pulse_ss = median( data_pulse ); % steady state after pulse
    
    data_depol = data_sweep(pulse_start_index:depol_stop_index);
    [depol_max,max_index] = max( data_depol ); % peak depol sag
    sag_depol = depol_max - pulse_ss;
    
    if plotnum == 1
        nSamples = size(data_sweep,1);
        dt = 1/Fs;
        time = (1:nSamples)*dt*1000;
        time_depol_peak = time( pulse_start_index + max_index );

        figure;
        plot(time, data_sweep,'-k','Linewidth', 2)
        hold on
        plot([time(1) time(end)],[pulse_ss pulse_ss],'--b')
        plot([time(1) time(end)],[depol_max depol_max],'--b')
        plot([time_depol_peak time_depol_peak],[depol_max pulse_ss],'--b')
        plot(time_depol_peak,depol_max,'ob')
        xlabel('Time (ms)')
        ylabel('Membrane Voltage (mV)')
        title( sprintf('Depol. SAG = %.1f mV',sag_depol) )
    end

end
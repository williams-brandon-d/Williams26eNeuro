function ISI_ratio = calculateISI(data_sweep,spike_indices,dt,plotnum)
        
    ISI_1 = (spike_indices(2) - spike_indices(1))*dt*1000;
    ISI_2 = (spike_indices(3) - spike_indices(2))*dt*1000;
    ISI_ratio = ISI_1/ISI_2;
    
    if plotnum == 1
        xlim_space = 5; % ms
        spike_line_min = -100; % mV
        spike_line_max = 100; % mV
        ISI_line_voltage = 0;
        ISI_color = {'b','g'};
        
        nSamples = size(data_sweep,1);
        time = (1:nSamples)*dt*1000;
        
        xlim_start = time(spike_indices(1)) - xlim_space;
        xlim_stop = time(spike_indices(3)) + xlim_space;
        xlimits = [xlim_start xlim_stop];

        figure;
        plot(time,data_sweep,'-k') % plot data
        hold on
        plot(time(spike_indices),data_sweep(spike_indices),'or') % plot spike peaks
        for iSpike = 1:3
            spike_index = spike_indices(iSpike);
            time_spike_index = time(spike_index);
            plot([time_spike_index time_spike_index],[spike_line_min spike_line_max],'--r') % plot spike time lines
        end
        for iISI = 1:2
            time_spike1 = time(spike_indices(iISI));
            time_spike2 = time(spike_indices(iISI+1));
            plot([time_spike1 time_spike2],[ISI_line_voltage ISI_line_voltage],'Linestyle','--','Color',ISI_color{iISI})
        end
        xlim(xlimits)
        xlabel('Time (ms)')
        ylabel('Membrane Voltage (mV)')
        title( sprintf( 'ISI 1/2 = %.2f', ISI_ratio) )
    end
    
end
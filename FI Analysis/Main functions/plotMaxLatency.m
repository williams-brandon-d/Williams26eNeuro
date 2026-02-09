function plotMaxLatency(data_sweep,spike_indices,pulse_start_index,dt)

    xlimits = [500 2500]; % Time (ms) for viewing only not analysis
    min_line_voltage = -100;
    max_line_voltage = 100;
    connect_line_voltage = 0;
    
    nSamples = size(data_sweep,1);
    time = (0:nSamples-1)*dt*1000; % time in msec
    time = time';
    
    first_spike_time_ms = time(spike_indices(1));
    pulse_start_time_ms = time(pulse_start_index);
    max_latency = first_spike_time_ms - pulse_start_time_ms;

    figure;
    plot(time,data_sweep,'-k')
    hold on
    if ~isempty(xlimits); xlim(xlimits); end
    plot([pulse_start_time_ms pulse_start_time_ms],[min_line_voltage max_line_voltage],'--r')
    plot([first_spike_time_ms first_spike_time_ms],[min_line_voltage max_line_voltage],'--r')
    plot([pulse_start_time_ms first_spike_time_ms],[connect_line_voltage connect_line_voltage],'--b')
    xlabel('Time (ms)')
    ylabel('Membrane Voltage (mV)')
    title( sprintf('Max Latency = %.2f ms',max_latency) )

end

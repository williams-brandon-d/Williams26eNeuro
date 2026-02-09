function [peak,auc] = getAmplitude(data,pulse_indices,dt,units,plotnum)
% write function for calculating peak amplitude and area under curve
% dt in sec

mean_data = mean(data(pulse_indices,:),2); % data is restricted to pulse stim indices 

mean_data = mean_data - mean_data(1); % first data point is baseline

peak = range(mean_data);

% instead of range use first peak height from baseline
% for inhibition peak is positive, excitation is negative

auc = trapz(mean_data)*dt*1000; % units*msec

if plotnum
    nSamples = numel(mean_data);
    time = (0:nSamples-1)*dt*1000; % time in msec
    plotAmplitude(mean_data,time,peak,auc,units);
end

function fig = plotAmplitude(data,time,peak,auc,units)
    fontsize = 15;
    fig = figure;
    plot(time,data,'k');
    title(sprintf('Peak = %d %s, AUC = %d %s*ms',round(peak),units,round(auc),units));
    xlabel('Time (ms)')
    ylabel(units)
    ax = gca;
    ax.XAxis.FontSize = fontsize;
    ax.XAxis.FontWeight = 'bold';
    ax.YAxis.FontSize = fontsize;
    ax.YAxis.FontWeight = 'bold';   
end


end
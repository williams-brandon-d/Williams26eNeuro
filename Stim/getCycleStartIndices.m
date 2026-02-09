function cycle_start_index = getCycleStartIndices(stim_raw, cycle_length, time, stim_plotnum )

    dt  = time(2) - time(1);
    Fs = 1/dt;
    stim_freq = round(Fs / cycle_length); % Hz
    
    lpf_stim = filterStimData(stim_raw,Fs,stim_freq);
    norm_stim = normalizeStimData(lpf_stim);  

%     % zero out pulse before theta stim
%     zeroEndTime = 0.5; % sec
%     maskTime = time < zeroEndTime;
%     norm_stim(maskTime) = min(norm_stim)*zeros(sum(maskTime),1);
%     % figure; plot(norm_stim)

    prom = (max(norm_stim)-min(norm_stim))/2;
    dist = 1000; % samples

    [~, min_locs] = findpeaks(-1*norm_stim,'MinPeakProminence',prom,'MinPeakDistance',dist);
    [~, max_locs] = findpeaks(norm_stim,'MinPeakProminence',prom,'MinPeakDistance',dist);

    stim_start_index = min_locs(1) - cycle_length;

    first_cycle_end = min_locs(1);
    first_cycle_start = first_cycle_end - cycle_length;
    last_cycle_end = min_locs(end) + cycle_length;
    nStimCycles = floor((last_cycle_end - first_cycle_start)/cycle_length);
    last_cycle_start = first_cycle_start + (nStimCycles-1)*cycle_length;

    cycle_start_index = first_cycle_start:cycle_length:last_cycle_start;

    if stim_plotnum == 1
        plotType = 2;
        min_peaks = norm_stim(min_locs);
        max_peaks = norm_stim(max_locs);
        min_time_peaks = time(min_locs);
        max_time_peaks = time(max_locs);
        switch plotType
            case 1
                figure; 
                plot(time,norm_stim,min_time_peaks,min_peaks,'ok',...
                    max_time_peaks,max_peaks,'or')
            case 2
                time_start = time(stim_start_index);
                stim_start = norm_stim(stim_start_index);
                time_windows_start = time(cycle_start_index);
                stim_windows_start = norm_stim(cycle_start_index);
                figure;
                plot(time,norm_stim,min_time_peaks,min_peaks,'ok',...
                    max_time_peaks,max_peaks,'or',...
                    time_start,stim_start,'xg',...
                    time_windows_start,stim_windows_start,'xm')
        end
    end
end
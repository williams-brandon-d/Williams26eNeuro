function lpf_stim = filterStimData(mean_stim_data,Fs,stim_freq)

    stim_lpf_fc = stim_freq*10; % Hz cutoff freq for filtering theta stim
    stim_lpf_order = 2;

    [b,a] = butter(stim_lpf_order,stim_lpf_fc/(Fs/2),'low'); % 
    lpf_stim = filtfilt(b,a,mean_stim_data);
end
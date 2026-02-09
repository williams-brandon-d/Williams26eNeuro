function norm_stim = normalizeStimData(stim_data)
    min_stim = min(stim_data);
    max_stim = max(stim_data);
    norm_stim = (stim_data - min_stim)/(max_stim - min_stim); % norm values from 0 to 1
end
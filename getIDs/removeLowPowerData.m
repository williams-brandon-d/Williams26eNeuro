function removeFlag = removeLowPowerData(file,dataSet,protocol)

fieldName = 'cell';

% check if fields exist in file
% if not, calculate

if isfield(file.(fieldName),'CWTmaxValues')
    phase = file.(fieldName).CWTmaxValues(1);
    freq = file.(fieldName).CWTmaxValues(2);
    power = file.(fieldName).CWTmaxValues(3);
else % calculate CWTmaxValues depending on protocol
    switch protocol
        case 'theta'
            
        case 'pulse'
            wavelet analysis parameters
            wname = 'amor'; % 'morse' (default), 'amor', 'bump'
            VoicesPerOctave = 32; % number of scales per octave
            flimits = [0 300]; % frequency limits for wavelet analysis
            [nSamples,~] = size(file.(fieldName).raw_data);
            setup filterbank for cwt 
            fb = cwtfilterbank('Wavelet',wname,'SignalLength',nSamples,...
            'FrequencyLimits',flimits,'SamplingFrequency',file.Fs,'VoicesPerOctave',VoicesPerOctave);
            gamma_data_noArtifacts = file.(fieldName).gamma_data(:,file.(fieldName).trials_noArtifacts);
            [x,y,z] = pulseCWT(gamma_data_noArtifacts,fb,file.time,file.CWTanalysisType);
            
            pulse_indices = file.pulse_start_index:file.pulse_stop_index;
            file.(fieldName).CWTmaxValues = getMaxValues(x(:,pulse_indices),y(:,pulse_indices),z(:,pulse_indices)); % max values during pulse
            phase = file.(fieldName).CWTmaxValues(1);
            freq = file.(fieldName).CWTmaxValues(2);
            power = file.(fieldName).CWTmaxValues(3);
    end
    save([file.saveFolder filesep 'data.mat'],'file','-mat','-nocompression'); % save updated file struct
end

if isfield(file.(fieldName),'CWTsnr')
    snr = file.(fieldName).CWTsnr;
else
    switch protocol
        case 'theta'
            
        case 'pulse'
        % wavelet analysis parameters
        wname = 'amor'; % 'morse' (default), 'amor', 'bump'
        VoicesPerOctave = 32; % number of scales per octave
        flimits = [0 300]; % frequency limits for wavelet analysis
        freq_threshold = 50; % zero frequencies below threshold for scalograms
        [nSamples,~] = size(file.(fieldName).raw_data);
        % setup filterbank for cwt 
        fb = cwtfilterbank('Wavelet',wname,'SignalLength',nSamples,...
        'FrequencyLimits',flimits,'SamplingFrequency',file.Fs,'VoicesPerOctave',VoicesPerOctave);
        gamma_data_noArtifacts = file.(fieldName).gamma_data(:,file.(fieldName).trials_noArtifacts);
        [x,y,z] = pulseCWT(gamma_data_noArtifacts,fb,file.time,file.CWTanalysisType);
        
        pulse_indices = file.pulse_start_index:file.pulse_stop_index;
        file.(fieldName).CWTmaxValues = getMaxValues(x(:,pulse_indices),y(:,pulse_indices),z(:,pulse_indices)); % max values during pulse
        if ~isnan(file.(fieldName).CWTmaxValues(1))
            file.(fieldName).CWTpeakStats = getXmaxPeakStats(x,y,z,file.(fieldName).CWTmaxValues(1),flimits,freq_threshold);
            file.(fieldName).CWTsnr = file.(fieldName).CWTmaxValues(3)/mean(z(:)); % peak power / average power
            snr = file.(fieldName).CWTsnr;
        else
            removeFlag = 1; % failed
            return;
        end
    end
    save([file.saveFolder filesep 'data.mat'],'file','-mat','-nocompression'); % save updated file struct
end

if isfield(file.(fieldName),'CWTpeakStats')
    width = file.(fieldName).CWTpeakStats.bandwidth;
else
    save([file.saveFolder filesep 'data.mat'],'file','-mat','-nocompression'); % save updated file struct
end

if ~isfield(file.(fieldName),'meanSpectralNoise')
    [file.(fieldName).meanSpectralNoise,file.(fieldName).stdSpectralNoise] = estimateSpectralNoise(file,fieldName,0); % estimate spectral noise
    save([file.saveFolder filesep 'data.mat'],'file','-mat','-nocompression'); % save updated file struct
end

% std_factor = 3; % peak data power must be 3 std above the mean spectral noise
% power_threshold = file.(fieldName).meanSpectralNoise + std_factor*file.(fieldName).stdSpectralNoise; % pA^2

switch dataSet
    case 'PV Transgenic'
        snr_threshold = 5;
        power_threshold = 10;
        width_threshold = 110;
    case 'Thy1'
        snr_threshold = 5;
        power_threshold = 20; % pA^2 (20)
        width_threshold = 100; % (Hz) bandwidth
    case 'Camk2'
        snr_threshold = 5;
        power_threshold = 20; % pA^2 (20)
        width_threshold = 100; % (Hz) bandwidth        
end

    freq_threshold = 55;

if strcmp(protocol,'pulse')
    width_threshold = inf;
    % phase = time in ms - must be within pulse indices
    pulse_indices = file.pulse_start_index:file.pulse_stop_index;
    pulse_time = file.time(pulse_indices);

    % conditions
    removeFlag = (freq > freq_threshold) & ...
                        (power > power_threshold) & ...
                        (phase > pulse_time(1)) & ...
                        (phase < pulse_time(end)) & ...
                        (snr > snr_threshold) & ...
                        (width < width_threshold);
    
    removeFlag = ~removeFlag; % failed at least one condition
%     removeFlag = 0; % dont remove
else
    phase_threshold = 3;
    
    % conditions
    removeFlag = (freq > freq_threshold) & ...
                        (power > power_threshold) & ...
                        (phase > -1*phase_threshold) & ...
                        (phase < phase_threshold) & ...
                        (snr > snr_threshold) & ...
                        (width < width_threshold);
    
    removeFlag = ~removeFlag; % failed at least one condition
end

end
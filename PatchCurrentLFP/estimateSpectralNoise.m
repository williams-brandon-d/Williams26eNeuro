function [meanNoise,stdNoise] = estimateSpectralNoise(file,fieldName,plotnum)
    
    % estimate scalogram noise during stim off time for ephys data
    
    wname = 'amor'; % 'morse' (default), 'amor', 'bump'
    VoicesPerOctave = 32; % number of scales per octave
    flimits = [0 300]; % frequency limits for wavelet analysis
    
    noise_length = 0.3; % seconds
    noise_start = 1; % start of file
    noise_stop = noise_length*file.Fs; % length of noise window
    
    noise_data = file.(fieldName).gamma_data(noise_start:noise_stop); % gamma filtered data
        
    sigLength = numel(noise_data);
    
    fb = cwtfilterbank('Wavelet',wname,'SignalLength',sigLength,...
    'FrequencyLimits',flimits,'SamplingFrequency',file.Fs,'VoicesPerOctave',VoicesPerOctave);    
    
    [wt,freq] = cwt(noise_data,'FilterBank',fb);
    z = abs( wt ).^2;

%     meanSpectralNoise = mean(z(:));
%     stdSpectralNoise = std(z(:));
    meanSpectralNoise = mean(z,2); % average spectrum across time (freq x 1)
    stdSpectralNoise = std(z,0,2); % average spectrum across time (freq x 1)

    maxDataFreq = file.(fieldName).CWTmaxValues(2);
    
    meanNoise = meanSpectralNoise(freq == maxDataFreq); % find mean of noise spectral power at max data frequency
    stdNoise = stdSpectralNoise(freq == maxDataFreq); % find std of noise spectral power at max data frequency

    if plotnum
        noise_time = (0:sigLength-1)*file.dt*1000; % ms
        [x,y] = meshgrid(noise_time,freq);
        figure; plot(noise_time,noise_data);
        plotScalogram(x,y,z,'Noise Spectrum',[],file.(fieldName).data_units,'noise');
    end
    
end
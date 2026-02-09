function [x,y,z,CWTcycleValues,coi] = thetaCWT(file,fieldName,cycle_start_index_noArtifacts,norm)

    % L2 normalization - signal energy is still not preserved due to numerical computation, but the result is "parseval-like"
    % wavelet magnitudes are not equal to signal amplitudes - scaling factor is dependent on frequency
    % L1 norm - signal amplitudes are equal to wavelet magnitudes
    % useful for interpretting scalograms
    if nargin < 4
        norm = 'L1'; % default matlab uses L1 wavelet normalization - signal energy is not preserved
    end

    % make scalogram for theta stim ephys data

    % 'amor' - the analytical morlet wavelet
    % center frequency: f_0 = 6 / (2*pi) Hz or w_0 = 6 rad/sec;
    % bandwidth: 2 Hz

    wname = 'amor'; % 'morse' (default), 'amor', 'bump'
    VoicesPerOctave = 32; % number of scales per octave
    flimits = [0 300]; % frequency limits for wavelet analysis
    
    delta_phase = 2*pi/file.cycle_length;
    cycle_phase = -pi:delta_phase:pi;

    N = numel(cycle_phase);

    fb = cwtfilterbank('Wavelet',wname,'SignalLength',N,...
    'FrequencyLimits',flimits,'SamplingFrequency',file.Fs,'VoicesPerOctave',VoicesPerOctave);

    % can compute frequencies and coi based on fb parameters
    freq = fb.centerFrequencies;
    [x,y] = meshgrid(cycle_phase,freq);

%     z = zeros(size(x)); % intialize scalogram

    % compute coi
%     Fc = 6 / (2*pi);
%     dt = 1/file.Fs;
%     time = (0:(N-1))*dt;
%     T_edge = min(time-time(1),time(end)-time); % τ[n]=min(t[n]−t[0],t[N−1]−t[n])
%     f_min = (sqrt(2) * Fc) ./ T_edge; % (1 x N) - shifted to the right by 1 compared to coi
% predtimes = (sqrt(2) * Fc) ./ freq;
    

    if strcmpi(norm,'L2')
        % f_0 = 6 / (2*pi) : center frequency of morlet wavelet
        % scales = (f_0 * file.Fs) ./ freq;
        scales = reshape(fb.scales,[],1); % (nFreq x 1)
    end

    
    % only use cycles without artifacts
    nCycles = numel(cycle_start_index_noArtifacts);
    
    cycles = 1:nCycles;

    CWTcycleValues = zeros(nCycles,3);
    
    for icycle = 1:numel(cycles)
           cycle = cycles(icycle);
           cycle_start = cycle_start_index_noArtifacts(cycle);
           cycle_stop = cycle_start + file.cycle_length; 
           
           window_data = file.(fieldName).gamma_data(cycle_start:cycle_stop);
    
%            figure;
%            subplot(2,1,1); plot(window_data)
%            subplot(2,1,2); plot(file.stim.raw_data(cycle_start:cycle_stop))

           if icycle == 1
               [wt,~,coi] = cwt(window_data,'FilterBank',fb);
               if strcmpi(norm,'L2')
                    wt = sqrt(scales).*wt;
               end
               z_cycle = abs( wt ).^2;
               z = z_cycle;
           else
               [wt,~,~] = cwt(window_data,'FilterBank',fb);
               if strcmpi(norm,'L2')
                    wt = sqrt(scales).*wt;
               end
               z_cycle = abs( wt ).^2;
               z = z + z_cycle;
           end

        CWTcycleValues(icycle,:) = getMaxValues(x,y,z_cycle);

    end
    z = z / nCycles;
        
end
function file = analyzePulseCurrentsLFP(fullfilename,printnum,savenum,saveFolder,info)

% analyze inhibition and excitation currents (plotPulseData,pulseCWT,pulsePSD,peak amplitude,sum (integrate) current during stim)
% analyze LFP (plotPulseData,pulseCWT,pulsePSD)

file.name = fullfilename;

file.info = info;

file.stimType = 'pulse';

file.CWTanalysisType = 'average scalogram'; % 'average scalogram' or 'average data'

% choose data to analyse
file.trials = 'all'; % choose trials for to average scalogram over

% filter params
file.gamma_bandpass = [50 250]; % Hz
file.filter_order = 4; % theta filter must be max order 4 - gamma can be higher

% wavelet analysis parameters
wname = 'amor'; % 'morse' (default), 'amor', 'bump'
VoicesPerOctave = 32; % number of scales per octave
flimits = [0 300]; % frequency limits for wavelet analysis
freq_threshold = 0; % zero frequencies below threshold for scalograms

% PSD parameters
% PSDtype = 'mtspec'; % 'mtspec' or 'pwelch'

% detect spike parameters
file.cell.range_threshold = 3000; % (pA) window data ranges > threshold are detected as spikes
file.lfp.range_threshold = 1000; % (uV) window data ranges > threshold are detected as spikes

% load data from .abf file
[data,si,file_info] = abfload(file.name,'start',0,'stop','e');

backslash_index = strfind(file_info.protocolName,'\');
file.protocol_name = file_info.protocolName(backslash_index(end)+1:end-4);

if size(data,2) == 3
    file.cell.raw_data = squeeze(data(:,1,:));
    file.lfp.raw_data = squeeze(data(:,2,:))*1000; % mV to uV
    file.stim.raw_data = squeeze(data(:,3,:));   
    file.cell.data_units = char(file_info.recChUnits(1));
    file.lfp.data_units = 'uV'; % lfp data is converted to uV instead of mV
    file.dataChannels = {'cell','lfp'};
else
    file.cell.raw_data = squeeze(data(:,1,:));
    file.stim.raw_data = squeeze(data(:,2,:));   
    file.cell.data_units = char(file_info.recChUnits(1));
    file.dataChannels = {'cell'};
end

file.nDataChannels = numel(file.dataChannels);

file.dt = si*(1e-6); % sampling interval (seconds)
file.Fs = 1/file.dt; % sampling frequency (Hz)

% get experimental params
file.led_input = file_info.DACEpoch.fEpochInitLevel(end); % mV
file.pulse_length = file_info.DACEpoch.lEpochPulsePeriod(end);

if isfield(file_info, 'comment')
    file.comment = file_info.comment;
else 
    file.comment = 'no comment found';
end

if printnum == 1
    fprintf('Protocol: %s\n',file.protocol_name);
    fprintf('LED Input = %g mV\n',file.led_input);
    fprintf('%s\n',file.comment);
end

[nSamples,~,ntrials] = size(data);

if strcmp(file.trials,'all')
   file.trials = 1:ntrials;
end
file.nTrials = numel(file.trials);

% construct time axis for plotting
time = (0:nSamples-1)*file.dt*1000; % time in msec
file.time = time'; % column vector

file.stim.norm = normalizeStimData(file.stim.raw_data); % normalize stim data from 0 to 1

% find stim start time
stim_threshold = 0.5;
mask1 = file.stim.norm(1:end-1) < stim_threshold;
mask2 = file.stim.norm(2:end) > stim_threshold;
file.pulse_start_index = find(mask1 & mask2);
file.pulse_start_time = file.time(mask1 & mask2); 

% find stim stop time
mask1 = file.stim.norm(1:end-1) > stim_threshold;
mask2 = file.stim.norm(2:end) < stim_threshold;
file.pulse_stop_index = find(mask1 & mask2);
file.pulse_stop_time = file.time(mask1 & mask2);

% setup filterbank for cwt 
fb = cwtfilterbank('Wavelet',wname,'SignalLength',nSamples,...
'FrequencyLimits',flimits,'SamplingFrequency',file.Fs,'VoicesPerOctave',VoicesPerOctave);

file.saveFolder = sprintf('%s - %g mV',saveFolder,file.led_input);
if ~exist(file.saveFolder, 'dir')
   mkdir(file.saveFolder)
end
    
% alternative add theta and gamma components to make threshold instead of
% using range threshold
file = detectArtifacts(file); % skip trials with artifacts

for i = 1:file.nDataChannels
    fieldName = file.dataChannels{i};

    % use trials without artifacts
    file.(fieldName).trials_noArtifacts = file.trials(~file.(fieldName).spikes); % assumes 1 trial of data
    file.(fieldName).nTrials_noArtifacts = numel(file.(fieldName).trials_noArtifacts);

    raw_data_noArtifacts = file.(fieldName).raw_data(:,file.(fieldName).trials_noArtifacts);

    %  baseline subtract instead of mean subtract - median during time before stim?
%     cell_data = file.(fieldName).raw_data - mean(file.(fieldName).raw_data,1); % mean subtract each trial

%     file.(fieldName).gamma_data = filterData(raw_data_noArtifacts,file.Fs,file.gamma_bandpass,file.filter_order);
    file.(fieldName).gamma_data = filterData(file.(fieldName).raw_data,file.Fs,file.gamma_bandpass,file.filter_order);

    gamma_data_noArtifacts = file.(fieldName).gamma_data(:,file.(fieldName).trials_noArtifacts);

    [x,y,z] = pulseCWT(gamma_data_noArtifacts,fb,file.time,file.CWTanalysisType);
    
    pulse_indices = file.pulse_start_index:file.pulse_stop_index;
    file.(fieldName).CWTmaxValues = getMaxValues(x(:,pulse_indices),y(:,pulse_indices),z(:,pulse_indices)); % max values during pulse
    if ~isnan(file.(fieldName).CWTmaxValues(1))
        file.(fieldName).CWTpeakStats = getXmaxPeakStats(x,y,z,file.(fieldName).CWTmaxValues(1),flimits,freq_threshold);
        file.(fieldName).CWTsnr = file.(fieldName).CWTmaxValues(3)/mean(z(:)); % peak power / average power
    end

    % plot all data trials and artifacts 
    figData = plotPulseData(file,fieldName);
    sgtitle(figData,sprintf( '%s - Pulse: 100 ms - Peak Gamma: %d Hz',file.name,round(file.(fieldName).CWTmaxValues(2)) ),'FontWeight','bold','Interpreter','none');
    
    % plot scalogram from data without artifacts 
    figCWT = plotScalogram(x,y,z,"",file.(fieldName).CWTmaxValues,file.(fieldName).data_units,file.stimType);
    
    % compute and plot PSD from data without artifacts
    [file.(fieldName).PSD,figPSD] = plotPulsePSD(raw_data_noArtifacts,file.Fs,file.gamma_bandpass,pulse_indices);

    if strcmp(fieldName,'cell')
        [file.(fieldName).amplitude_peak,file.(fieldName).amplitude_sum] = getAmplitude(file.(fieldName).raw_data,pulse_indices,file.dt,file.(fieldName).data_units,0);
        format longG
        fprintf('%s Range = %d %s\n',fieldName,round(file.(fieldName).amplitude_peak),file.(fieldName).data_units);
        fprintf('%s Area Under Curve = %d %s*ms\n',fieldName,round(file.(fieldName).amplitude_sum),file.(fieldName).data_units);
    end

    if printnum == 1
        fprintf('%s Gamma Frequency: %d Hz\n',fieldName,round(file.(fieldName).CWTmaxValues(2)));
    end
    
    if savenum == 1
%         saveas(figData,[file.saveFolder filesep fieldName '.svg']);
%         saveas(figCWT,[file.saveFolder filesep fieldName ' scalogram.svg']);
%         saveas(figPSD,[file.saveFolder filesep fieldName ' Raw Data Trial PSD.svg']);
        print(figData,'-vector','-dsvg',[file.saveFolder filesep fieldName '.svg']);
        print(figCWT,'-vector','-dsvg',[file.saveFolder filesep fieldName ' scalogram.svg']);
        print(figPSD,'-vector','-dsvg',[file.saveFolder filesep fieldName ' Raw Data Trial PSD.svg']);
    end
    
end

if file.nDataChannels == 2
    [file,figXcorr] = pulseLFPxcorr(file); % cycle current-lfp cross corr
    if savenum == 1
        saveas(figXcorr,[file.saveFolder filesep fieldName ' Xcorr.svg']);
    end
end


end

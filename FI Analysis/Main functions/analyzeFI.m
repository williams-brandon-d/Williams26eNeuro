function file = analyzeFI(fullfilename,printnum,savenum,saveFolder,info)

file.name = fullfilename;
file.info = info;
file.saveFolder = [saveFolder ' FI'];

if ~exist(file.saveFolder, 'dir')
   mkdir(file.saveFolder)
end

cellType = info.cell_type; % cell type hypothesis

% plot flags | 0 = no plot, 1 = plot data
sag_plotnum = 0;
ISI_plotnum = 0;
max_depol_plotnum = 0;
max_latency_plotnum = 0; 
max_dap_plotnum = 0; % just max sweep
data_plotnum = 0;
tau_plotnum = 0;
dap_plotnum = 0; % each sweep
depol_sag_plotnum = 0; % each sweep

% spike detection params
spike_threshold = 0; % mV
spike_dist_ms = 1; % ms

% low pass data filter params
data_lpf_fc = 2000; % Hz cutoff freq for window data
data_lpf_order = 8;

% plot data params
plot_sweeps = [];
plotTitle = ''; % if empty, auto generates title


% read abf file
[data,si,file_info] = abfload(fullfilename,'start',0,'stop','e');

% parse file info
backslash_index = strfind(file_info.protocolName,'\');
file.protocol_name = file_info.protocolName(backslash_index(end)+1:end-4);

% protocol current pulse start and stop time
switch file.protocol_name
    case 'FI_30currentsteps_FRF'
        % current = -100:25:525 | 26 sweeps
        % 2 seconds off, 1 second pulse, 2 seconds off
        pulse_start_sec = 1.975; % pulse starts before 2 sec
    otherwise
        % current_sweeps = -200:25:525; % pA
        % 1 seconds off, 1 second pulse, 1 seconds off
        pulse_start_sec = 1; % pulse starts at 1 sec
end

pulse_stop_sec = pulse_start_sec + 1;

file.data_units = char(file_info.recChUnits(1));

if isfield(file_info, 'comment')
    file.comment = file_info.comment;
else 
    file.comment = 'no comment found';
end

% get current sweep values
file.current_min = file_info.DACEpoch.fEpochInitLevel(2); % data units
file.current_delta = file_info.DACEpoch.fEpochLevelInc(2); % data units

% not accurate
% pulse_start_index = file_info.DACEpoch.lEpochInitDuration(1) + 1; % samples
% pulse_stop_index = pulse_start_index + file_info.DACEpoch.lEpochInitDuration(2) - 1; % samples

% dt and Fs
dt = si*(1e-6); % sampling interval (seconds)
Fs = 1/dt; % sampling frequency (Hz)

% analyze data
data = squeeze(data);
nSweeps = size(data,2);

file.current_max = (nSweeps-1)*file.current_delta + file.current_min;

% filter data
[b_lpf,a_lpf] = butter(data_lpf_order, data_lpf_fc/(Fs/2), 'low');
data = filtfilt(b_lpf, a_lpf, data); 

pulse_start_index = uint64(Fs*pulse_start_sec + 1);
pulse_stop_index = uint64(Fs*pulse_stop_sec);

% calculate hyperpolarizing sag from most hyperpolarizing current pulse
% calculate depolarizing sag from the rebound of the same pulse
[sag_hyper,sag_depol] = calculateSAG(data,pulse_start_index,pulse_stop_index,dt,sag_plotnum);

spike_dist_sec = spike_dist_ms/1000;
spike_point_distance = Fs*spike_dist_sec;

sag_depol_thresh = cell(nSweeps,1);
latency = cell(nSweeps,1);
average_dap = cell(nSweeps,1);
ISI_ratio = cell(nSweeps,1);
perithreshold = cell(nSweeps,1);

dap_count = 0;
% perithreshold_sweep = [];

current_sweeps = file.current_min:file.current_delta:file.current_max;

% calculate resting membrane potential
restMask = current_sweeps == 0; 
file.RMP = median(data(:,restMask));

% could run spike detection on all sweeps first 

for iSweep = 1:nSweeps
    data_sweep = data(:, iSweep);
    
    if (max(data_sweep) > spike_threshold) % prominence and width instead of threshold?
        % detect spikes
        [~, spike_indices] = findpeaks(data_sweep,'MinPeakHeight',spike_threshold,'MinPeakDistance',spike_point_distance);
        nSpikes = length(spike_indices);
        
        % calculate DAP
        if current_sweeps(iSweep) < 300 && dap_count < 4
            [average_dap{iSweep},~] = calculateDAP(data_sweep,spike_indices,Fs,pulse_start_index,pulse_stop_index,dap_plotnum);
            if ~isempty(average_dap{iSweep}); dap_count = dap_count + 1; end
        end
        
        % calculate latency to first spike
        if (spike_indices(1) < pulse_stop_index) && (spike_indices(1) > pulse_start_index)
            latency{iSweep} = (spike_indices(1)-double(pulse_start_index))*dt*1000;
        end
        
        % calculate ISI ratio
        if nSpikes >= 3
            ISI_ratio{iSweep} = calculateISI(data_sweep,spike_indices,Fs,ISI_plotnum);
        end 
        
    elseif current_sweeps(iSweep) > 0
        % calculate depolarizing sag
        sag_depol_thresh{iSweep} = calculateDepolSAG(data_sweep,Fs,pulse_start_index,pulse_stop_index,depol_sag_plotnum);
        
        % save perithreshold sweep index for time constant calculation
        if max(data_sweep(pulse_stop_index:end)) < spike_threshold % check for spikes after pulse
%             perithreshold_sweep = iSweep;
            % instead find sweep with largest change in potential during pulse vs after
            pulse_median = median(data_sweep(pulse_start_index:pulse_stop_index));
            pulse_after = median(data_sweep(pulse_stop_index:end));
            perithreshold{iSweep} = pulse_median - pulse_after;
        end
    end
end

% perithreshold sweep
peri_not_empty = find(~cellfun(@isempty,perithreshold));
perithreshold_mat = cell2mat(perithreshold);
[~,max_perithreshold_mat_index] = max(perithreshold_mat);
perithreshold_sweep = peri_not_empty(max_perithreshold_mat_index);

% time constant
if ~isempty(perithreshold_sweep)
    file.time_constant = calculateTimeConstant(data(:,perithreshold_sweep),Fs,pulse_start_index,pulse_stop_index,tau_plotnum,cellType);
else
    file.time_constant = NaN;
end

% max depolarizing sag potential
depol_not_empty = find(~cellfun(@isempty,sag_depol_thresh));
sag_depol_mat = cell2mat(sag_depol_thresh);
[sag_depol_threshold,max_depol_mat_index] = max(sag_depol_mat);
max_depol_index = depol_not_empty(max_depol_mat_index);
if max_depol_plotnum == 1 && ~isempty(max_depol_index)
    data_sweep = data(:,max_depol_index);
    calculateDepolSAG(data_sweep,Fs,pulse_start_index,pulse_stop_index,max_depol_plotnum);
end

% max latency 
latency_not_empty = find(~cellfun(@isempty,latency));
latency_mat = cell2mat(latency);
[max_latency,max_latency_mat_index] = max(latency_mat);
max_latency_index = latency_not_empty(max_latency_mat_index);
if max_latency_plotnum == 1
   data_sweep = data(:,max_latency_index);
   [~, spike_indices] = findpeaks(data_sweep,'MinPeakHeight',spike_threshold,'MinPeakDistance',spike_point_distance);
   plotMaxLatency(data_sweep,spike_indices,pulse_start_index,dt)
end

% max DAP
dap_not_empty = find(~cellfun(@isempty,average_dap));
dap_mat = cell2mat(average_dap);
[max_dap,max_dap_mat_index] = max(dap_mat);
if max_dap == 0
    max_dap_index = max_latency_index;
else
    max_dap_index = dap_not_empty(max_dap_mat_index);
end
if max_dap_plotnum == 1
    data_sweep = data(:,max_dap_index);
    [~, spike_indices] = findpeaks(data_sweep,'MinPeakHeight',spike_threshold,'MinPeakDistance',spike_point_distance);
    [~,~] = calculateDAP(data_sweep,spike_indices,Fs,pulse_start_index,pulse_stop_index,max_dap_plotnum);
end

% plot data example traces
if data_plotnum == 1 
    if isempty( plot_sweeps )
        plot_sweeps = 1:nSweeps;
    elseif strcmp(plot_sweeps,'all')
        plot_sweeps = 1:nSweeps;            
    elseif max(plot_sweeps) > nSweeps
        if mod(nSweeps,2) == 1
            plot_sweeps = 1:2:nSweeps;
        else
            plot_sweeps = 1:2:(nSweeps-1);
        end
    end
    
    plotFIdata(data,plot_sweeps,dt,plotTitle,info)
end 

ISI_array = cell2mat(ISI_ratio);

% print feature results
if printnum == 1
    fprintf('Hyper Sag = %.2f mV\n',sag_hyper)
    fprintf('Depol Sag = %.2f mV\n',sag_depol)
    fprintf('Threshold Depol Sag = %.2f mV\n',sag_depol_threshold)
    fprintf('First ISI 1/2 = %.2f\n',ISI_array(1))
    fprintf('Second ISI 1/2 = %.2f\n',ISI_array(2))
    fprintf('Latency to First Spike = %.2f ms\n',max_latency)
    fprintf('dAP = %.2f mV\n',max_dap)
end

% save data in file struct
file.sag_depol = sag_depol;
if sag_depol_threshold < 0
    sag_depol_threshold = 0;
end
file.sag_depol_thresh = sag_depol_threshold;
file.sag_hyper = sag_hyper;
file.dap = max_dap;
file.latency = max_latency;

%     if ISI_array(1) > 2
%         info(ID_index).ISI_ratio = ISI_array(2);
%     else
%         info(ID_index).ISI_ratio = ISI_array(1);
%     end

if numel(ISI_array) > 2 
    file.ISI_ratio = median(ISI_array(1:3));
else
    file.ISI_ratio = median(ISI_array(1:end));
end
%     fprintf('Median ISI 1/2 = %.2f\n',info(ID_index).ISI_ratio)

file.dt = dt;
file.Fs = Fs;
file.pulse_start_index = pulse_start_index;
file.pulse_stop_index = pulse_stop_index;



end
function plotLightvsCellPaired()

clear variables; close all; clc;
cd('C:\Users\brndn\OneDrive\Desktop\White Lab\Matlab Files\');
addpath(genpath('./')); % add all folders and subfolders in cd to path

savenum = 1; % 0 = don't save, 1 = save info

% plot light power vs gamma power/freq, average firing rate

% test nonlinear monotonic correlation - spearman??
% [RHO,PVAL] = corr(a',b','Type','Spearman');

% PV Transgenic dataset 'theta' protocol is mostly single light intensities
% try 'pulse' for PV Transgenic dataset
% problem with peak frequency analysis?

dataSets = {'PV Transgenic'};
cell_types = {'stellate','pyramidal','fast spiking'};
% cell_types = {'stellate','pyramidal'};

% dataTypes = {'power','freq','rate'};
dataTypes = {'power','freq'};

protocol = 'pulse'; % theta or pulse

% find ID parameters for analysis
params.locations = 'all';
params.cell_nums = 'all';
params.comments = {'','DNQX before','DNQX before 10 uM','GBZ before','DNQX_GABAzine_before', 'Gabazine before'...
                   'x','DNQX beforex','DNQX before 10 uMx','GBZ beforex','DNQX_GABAzine_beforex', 'Gabazine beforex'}; % 'GBZ before' 'DNQX_GABAzine_before', 'Gabazine before', 'DNQX before'

switch protocol
    case 'theta'
        params.protocols = {'theta','theta_2chan'};
    case 'pulse'
        params.protocols = {'pulse','pulse_2chan'};
end

if strcmp(dataSets,'all'); dataSets = {'Camk2','Thy1','PV Transgenic','PV Viral'}; end

nCellTypes = numel(cell_types);
nDataTypes = numel(dataTypes);

for iSet = 1:numel(dataSets)
dataSet = dataSets{iSet};

[info,~,data_path] = getInfo(dataSet);

switch dataSet
    case 'Thy1'
         saveFolder = 'C:\Users\brndn\Downloads\Thy1-ChR2\Raw Data\mEC\results\Summary';
%          experiments = {'inhibition'};
    case 'PV Transgenic'
         saveFolder = 'C:\Users\brndn\Downloads\PV-ChR2 Transgenic\Summary';
    case 'PV Viral'
         saveFolder = 'C:\Users\brndn\Downloads\PV-ChR2\Summary';
    case 'Camk2'
         saveFolder = 'C:\Users\brndn\Downloads\CaMK2-ChR2\Summary';
%          experiments = {'excitation'};
end

for iType = 1:nDataTypes
dataType = dataTypes{iType};

switch dataType
    case 'power'
        experiment = 'inhibition';
    case 'freq'
        experiment = 'inhibition';
    case 'rate'
        experiment = 'currentclamp';
end

params.experiments = {experiment};

for iCell = 1:nCellTypes
params.cell_types = cell_types(iCell);
cellType = cell_types{iCell};

IDs = getIDs(info,params);

IDs = removeIDs(IDs,info,dataSet); % skip bad recordings - get params for cells to skip

% validIDs = removeBadLFP(IDs,info,dataSet); % remove bad LFP recordings (CamK2 only)

if (isempty(IDs)); fprintf('No %s Files Found.',cell_types{iCell}); continue; end

nIDs = numel(IDs);

ID_index = zeros(nIDs,1);
for m = 1:nIDs
    ID_index(m) = find_index(info,'ID',IDs{m});
end

info_IDs = info(ID_index);

cellNums = {info_IDs.cell_num};
diffNums = unique(cellNums);
[~,lengthIdx] = sort(cellfun(@length,diffNums),'ascend');
uniqueCellNums = diffNums(lengthIdx);
nCellNums = numel(uniqueCellNums);


fig = figure;
hold on;

fontsize = 12;
fontweight = 'bold';
ax = gca;
ax.YAxis.FontSize = fontsize;
ax.XAxis.FontSize = fontsize;
ax.YAxis.FontWeight = fontweight;
ax.XAxis.FontWeight = fontweight;
axis square;

% change x and y axis limits for PV and pulse data
switch dataSet
    case 'PV Transgenic'
        switch protocol
            case 'theta'
                xlimits = [0 15]; % mW / mm^2
            case 'pulse'
                switch cellType
                    case 'fast spiking'
                        xlimits = [0 15]; % mW / mm^2
                    case 'pyramidal'
                        xlimits = [0 4];
                    case 'stellate'
                        xlimits = [0 4];
                end
        end
    otherwise
        xlimits = [0 25]; % mW / mm^2
end

xlim(xlimits);
xlabel('Light Intensity (mW/mm^{2})','FontSize',12,'FontWeight','bold');

switch dataType
    case 'power'
        switch dataSet
            case 'PV Transgenic'
                switch protocol
                    case 'theta'
                        ylimits = [0 4]; % log power
                    case 'pulse'
                        ylimits = [1 5]; % log power
                end
            otherwise
                ylimits = [1 5]; % log power
        end
        ystring = sprintf('%s Peak Gamma Power (Log pA^{2})',dataSet);
    case 'freq'
        switch dataSet
            case 'PV Transgenic'
                switch protocol
                    case 'theta'
                        ylimits = [50 200]; % Hz
                    case 'pulse'
                        ylimits = [50 250]; % Hz
                end
            otherwise
                ylimits = [50 200]; % Hz
        end
        ystring = sprintf('%s Peak Frequency (Hz)',dataSet);
    case 'rate'
        ylimits = [0 5]; % Hz
        ystring = sprintf('%s Firing Rate (Spikes / Theta Cycle)',dataSet);      
end

ylim(ylimits);
ylabel(ystring,'FontSize',12,'FontWeight','bold');

for iNum = 1:nCellNums
    cellnum = uniqueCellNums{iNum};
    cellnum_mask = ismember(cellNums, cellnum); % find info_IDs indices for cellnum
    cell_info = info_IDs(cellnum_mask);

    nCellIDs = numel(cell_info);

    % sort info by light power
    xArray = zeros(nCellIDs,1);
    yArray = zeros(nCellIDs,1);
    gammaPowers = zeros(nCellIDs,1);

    for iCellID = 1:nCellIDs
    
        p = cell_info(iCellID);

        ID = p.ID;

        filename = sprintf('%s.abf',ID);
        fprintf('Analyzing %s,File %d/%d\n',filename,iCellID,nCellIDs);

        if isempty(p.comments)
            comment = 'No Comment';
        else
            comment = p.comments;
        end
    
        commentFolder = sprintf('%sresults\\%s\\%s\\%s\\%s\\%s\\%s\\%s',data_path,p.location,p.cell_type,p.cell_num,p.experiment,p.protocol,comment);
    
        dataFolder = getIDFolder(commentFolder,ID);
    
        dataFilename = [dataFolder filesep 'data.mat'];
    
        load(dataFilename,'file');

        gammaPowers(iCellID) = log10(file.cell.CWTmaxValues(3));
    
        % gather relevant data 
        if ~isempty(file.led_input)
            focalArea = 0.237; % mm^2 40x objective 
            % convert led input to light intensity
            power = LEDvoltageToPower470nm(file.led_input); % mW
            intensity = power / focalArea; % mW / mm^2
%             inputArray(iID) = file.led_input;
            xArray(iCellID) = intensity;

%             fieldName = 'cell';
%             % wavelet analysis parameters
%             wname = 'amor'; % 'morse' (default), 'amor', 'bump'
%             VoicesPerOctave = 32; % number of scales per octave
%             flimits = [0 300]; % frequency limits for wavelet analysis
%             freq_threshold = 50; % zero frequencies below threshold for scalograms
%             [nSamples,~] = size(file.(fieldName).raw_data);
%             % setup filterbank for cwt 
%             fb = cwtfilterbank('Wavelet',wname,'SignalLength',nSamples,...
%             'FrequencyLimits',flimits,'SamplingFrequency',file.Fs,'VoicesPerOctave',VoicesPerOctave);
%             gamma_data_noArtifacts = file.(fieldName).gamma_data(:,file.(fieldName).trials_noArtifacts);
%             [x,y,z] = pulseCWT(gamma_data_noArtifacts,fb,file.time,file.CWTanalysisType);
%             
%             pulse_indices = file.pulse_start_index:file.pulse_stop_index;
%             file.(fieldName).CWTmaxValues = getMaxValues(x(:,pulse_indices),y(:,pulse_indices),z(:,pulse_indices)); % max values during pulse
%             file.(fieldName).CWTpeakStats = getXmaxPeakStats(x,y,z,file.(fieldName).CWTmaxValues(1),flimits,freq_threshold);
%             file.(fieldName).CWTsnr = file.(fieldName).CWTmaxValues(3)/mean(z(:)); % peak power / average power
%             save([file.saveFolder filesep 'data.mat'],'file','-mat','-nocompression'); % save updated file struct

            switch dataType 
                case 'power'
                    yArray(iCellID) = log10(file.cell.CWTmaxValues(3)); % log lfp peak scalogram power
                case 'freq'
                    % remove frequencies with low power etc
                    removeFlag = removeLowPowerData(file,dataSet,protocol);
                    if removeFlag
                        yArray(iCellID) = NaN;
                    else 
                        yArray(iCellID) = file.cell.CWTmaxValues(2); % lfp peak scalogram freq
                    end
%                     yArray(iCellID) = file.cell.CWTmaxValues(2); % lfp peak scalogram freq
                case 'rate'
                    % calculate average firing rate
    %                 file.FRs = FRs;
                    nSpikes = numel([file.theta_spike_phases{:}]); % nCycles x cycle phase array
                    nCycles = numel(file.cycles);
                    if ~isempty(nSpikes)
                        yArray(iCellID) = nSpikes / nCycles; % average spikes per theta cycle 
                    else
                        yArray(iCellID) = 0; % change to NaN to ignore non-firing cells
                    end
            end
        else
            xArray(iCellID) = NaN; % ignore if either value is empty
            yArray(iCellID) = NaN;
        end
    
        clearvars file;
    end

    % sort by light intensity
    [xArray,sortedIdx] = sort(xArray,'ascend');
    yArray = yArray(sortedIdx);

    % - remove duplicates of same light intensity - keep greater gamma power
    uniqueLight = unique(xArray); % get unique values
    if numel(uniqueLight) ~= numel(xArray)
        for iLight = 1:numel(uniqueLight)
            % find duplicates
            light = uniqueLight(iLight);
            indices = find(xArray == light);
            % duplicate gamma powers
            lightGammaPower = gammaPowers(indices);
            [~,sortedGammaPowerIdx] = sort(lightGammaPower,'descend'); % sorted
            % remove data besides highest power from xArray and yArray
            removeIndices = indices(sortedGammaPowerIdx(2:end));
            xArray(removeIndices) = [];
            yArray(removeIndices) = [];
        end
    end

    switch cellType
        case 'stellate'
            color = [1 0 0];
        case 'pyramidal'
            color = [0 1 0];
        case 'fast spiking'
            color = [0 0 1];
    end
    
    %     plot(xCell{iCell},yCell{iCell},'.','Color',color,'DisplayName',cellType)
    plot(xArray,yArray,'Color',color,'Marker','.','MarkerSize',10,'DisplayName',cellType);
    plot(xArray,yArray,'Color',color);

end 

if savenum
%     saveFolder = 'C:\Users\brndn\OneDrive\Desktop\White Lab\Carmen Grants';
    saveFilename = [saveFolder filesep sprintf('%s light power vs %s Paired %s %s.svg',dataSet,dataType,protocol,cellType)];
    print(fig,'-vector','-dsvg',saveFilename);
end

end

end

end

close all;

    function powers = LEDvoltageToPower470nm(voltages)
        % measured LED input voltages and light powers
        measuredVoltages =    [0.500 0.600 0.700 0.800 0.900 1.000 1.500 2.000 5.000 10.000]; % mV
        measuredPowers = [0.000 0.040 0.078 0.110 0.145 0.180 0.350 0.495 1.425 2.8450]; % mW

        % linear interpolation and extrapolation
        powers = interp1(measuredVoltages,measuredPowers,voltages,"linear",'extrap');

    end

end
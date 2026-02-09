function boxplotCellData(dataSets, comment_group1, savenum)
% independent group plots and comparisons

% dataTypes = {'power','frequency','phase'};
dataTypes = {'totalGammaPower','totalFastPower'};
% dataTypes = {'frequencyACF'};

% for total power
gammaRange = [60 140];
fastRange = [100 200];

wavelet_norm = 'L2'; % L1 (default) or L2

sig = 'symbol'; % 'symbol' or 'exact' | p-value significance indicator

params.comments = comment_group1;
comments = cell2mat(params.comments);

% find ID parameters for analysis
params.locations = 'all';
params.cell_nums = 'all';
params.protocols = {'theta','theta_2chan'};
protocol = params.protocols{1};

if strcmp(dataSets,'all'); dataSets = {'Camk2','Thy1','PV Transgenic','PV Viral'}; end

nDataTypes = length(dataTypes);

for iSet = 1:numel(dataSets)
dataSet = dataSets{iSet};

[info,~,data_path] = getInfo(dataSet);

switch dataSet
    case 'Thy1'
         saveFolder = 'C:\Users\brndn\Downloads\Thy1-ChR2\Raw Data\mEC\results\Summary';
         experiments = {'excitation','inhibition'};
%          experiments = {'inhibition'};
    case 'PV Transgenic'
         saveFolder = 'C:\Users\brndn\Downloads\PV-ChR2 Transgenic\Summary';
         experiments = {'inhibition'};
    case 'PV Viral'
         saveFolder = 'C:\Users\brndn\Downloads\PV-ChR2\Summary';
         experiments = {'inhibition'};
    case 'Camk2'
         saveFolder = 'C:\Users\brndn\Downloads\CaMK2-ChR2\Summary';
         experiments = {'excitation','inhibition'};
%          experiments = {'excitation'};
end

if ~exist(saveFolder, 'dir')
   mkdir(saveFolder)
end

nExperiments = length(experiments);

for iExp = 1:nExperiments
params.experiments = experiments(iExp);
experiment = experiments{iExp};

switch experiment
    case 'inhibition'
        cell_types = {'stellate','pyramidal','fast spiking'};
%         cell_types = {'stellate'};
    case 'excitation'
        cell_types = {'stellate','pyramidal','fast spiking'};
%         cell_types = {'fast spiking'};
end

nCellTypes = length(cell_types);

allFreqACF = cell(nCellTypes,1);
allGammaPower = cell(nCellTypes,1);
allFastPower = cell(nCellTypes,1);
allPower = cell(nCellTypes,1);
allFreq = cell(nCellTypes,1);
allPhase = cell(nCellTypes,1);

for iCell = 1:nCellTypes
    cellType = cell_types{iCell};
    params.cell_types = cell_types(iCell);

    IDs = getIDs(info,params);

    IDs = removeIDs(IDs,info,dataSet); % skip bad recordings - get params for cells to skip
    
    if (isempty(IDs)); fprintf('No %s %s Files Found.',params.cell_types,comments); continue; end
    
    nIDs = numel(IDs);

    freqACFArray = cell(nIDs,1);
    gammaPowerArray = cell(nIDs,1);
    fastPowerArray = cell(nIDs,1);
    powerArray = cell(nIDs,1);
    freqArray = cell(nIDs,1);
    phaseArray = cell(nIDs,1);

    for iID = 1:nIDs
    
        ID = IDs{iID};
        ID_index = find_index(info,'ID',ID);

        p = info(ID_index);
    
        if isempty(p.comments)
            comment = 'No Comment';
        else
            comment = p.comments;
        end

        filename = sprintf('%s.abf',ID);
        fprintf('%s %s %s Analyzing %s,File %d/%d\n',dataSet,cellType,comment,filename,iID,nIDs)
    
        commentFolder = sprintf('%sresults\\%s\\%s\\%s\\%s\\%s\\%s\\%s',data_path,p.location,p.cell_type,p.cell_num,p.experiment,p.protocol,comment);
    
        dataFolder = getIDFolder(commentFolder,ID);
    
        dataFilename = [dataFolder filesep 'data.mat'];
    
        load(dataFilename,'file');        

        % gather relevant data 
        % peak power values
        phase = file.cell.CWTmaxValues(1);
        freq = file.cell.CWTmaxValues(2);
        power = file.cell.CWTmaxValues(3);

        if any(strcmp(dataTypes,{'totalGammaPower','totalFastPower'}))
            % total power - compare L1 and L2 normalization
            fieldName = 'cell';
            [x,y,z,file.(fieldName).CWTcycleValues,coi] = thetaCWT(file,fieldName, file.cycle_start_index_noArtifacts,wavelet_norm);
            file.(fieldName).CWTmaxValues = getMaxValues(x,y,z);
            file.(fieldName).CWTpeakStats = getXmaxPeakStats(x,y,z,file.(fieldName).CWTmaxValues(1),[0 200],0);
            file.(fieldName).CWTsnr = file.(fieldName).CWTmaxValues(3)/mean(z(:)); % peak power / average power
    
            % only sum over regions in the cone of influence (COI) - ignores boundary effects of wavelet transform
            coi = reshape(coi,1,[]); % (1 x N) % coi represents the lowest valid frequency for each time point
            coiMask = y >= coi;   % results in [numFreq x N] logical mask
            z = z.*coiMask; % zero wavelet coefficients outside COI
    
            f = y(:,1); % frequency axis of scalogram (nFreqs x 1)
    
            % sum power over different frequency ranges
            gammaMask = (f >= gammaRange(1) & f <= gammaRange(2));
            totalGammaPower = sum(sum(z(gammaMask,:)));  % integrate over gamma frequencies and time
    
            fastMask = (f >= fastRange(1) & f <= fastRange(2));
            totalFastPower = sum(sum(z(fastMask,:)));  % integrate over gamma frequencies and time

            gammaPowerArray{iID} = totalGammaPower;
            fastPowerArray{iID} = totalFastPower;
        end

        % if data is frequency or phase - remove low power data 
        removeFlag = removeLowPowerData(file,dataSet,protocol);
        if removeFlag; freq = []; phase = []; end

        if any(strcmp(dataTypes,'frequencyACF'))
            % --- Estimate with autocorrelation
            plotnum = 0;
            cycles = 2:11;
            nCycles = numel(cycles);
            acfFreqs = zeros(nCycles,1);
            peakValues = acfFreqs;
            welchValues = acfFreqs;
            for iCycle = 1:nCycles
               cycle = cycles(iCycle);
               cycle_start = file.cycle_start_index_noArtifacts(cycle);
               cycle_stop = cycle_start + file.cycle_length; 
               
               window_data = file.cell.gamma_data(cycle_start:cycle_stop);
                [acfFreqs(iCycle), ~, infoACF,~] = acf_freq_estimate(window_data, file.Fs, plotnum, ...
                    'FreqRange', [50 200], ...        % Hz (set empty [] to disable)
                    'Bandpass',  [], ...        % Hz (set [] to skip)
                    'Detrend',   false, ...
                    'MaxLag',    0.02, ...            % seconds of ACF to compute
                    'MinPeakProminence', 0); % corr: first positive peak - first negative peak
                peakValues(iCycle) = infoACF.peakValue;
                welchValues(iCycle) = infoACF.fftpeak;
    %                         fprintf('ACF estimate: f0 = %.2f Hz, peak r = %.3f\n', acfFreqs(iCycle), peakValues(iCycle));
            end
            medACF = median(acfFreqs,'omitnan');
            meanR = mean(peakValues,'omitnan');
            medFFT = median(welchValues,'omitnan');
            fprintf('ACF estimate: f0 = %.2f Hz, peak r = %.3f, Welch estimate: f = %.2f Hz\n', medACF, meanR, medFFT);
    %                     data = medFFT;
            freqACFArray{iID} = medACF;
        end
        
        powerArray{iID} = power;
        freqArray{iID} = freq;
        phaseArray{iID} = phase;

        clearvars file;
    end

    allFreqACF{iCell,1} = cell2mat(freqACFArray); % save cell data
    allGammaPower{iCell,1} = cell2mat(gammaPowerArray); % save cell data
    allFastPower{iCell,1} = cell2mat(fastPowerArray); % save cell data
    allPower{iCell,1} = cell2mat(powerArray); % save cell data
    allFreq{iCell,1} = cell2mat(freqArray); % save cell data
    allPhase{iCell,1} = cell2mat(phaseArray); % save cell data
end


for iType = 1:nDataTypes
    dataType = dataTypes{iType};
    
    switch dataType
        case 'frequencyACF'
            allData = allFreqACF;
        case 'totalGammaPower'
            allData = allGammaPower;
        case 'totalFastPower'
            allData = allFastPower;
        case 'power'
            allData = allPower;
        case 'frequency'
            allData = allFreq;
        case 'phase'
            allData = allPhase;
    end

    saveFilename = [saveFolder filesep sprintf('%s %s %s %s stats.xlsx',dataSet,protocol,experiment,dataType)];
    stats = myMultipleIndependentGroupStats(allData,cell_types,comments,'linear',saveFilename,savenum);

    if nCellTypes == 1
        stats.kw = [];
    end

    stats.kw.sig = sig;
    
    % convert data to log for power plots
    if any(strcmp(dataType,{'power','totalFastPower','totalGammaPower'}))
        saveFilename = [saveFolder filesep sprintf('%s %s %s %s LOG stats.xlsx',dataSet,protocol,experiment,dataType)];
        LOGstats = myMultipleIndependentGroupStats(allData,cell_types,comments,'log',saveFilename,savenum);
        % log transform data?
         allData = cellfun(@log10,allData,'UniformOutput',false); % returns cell array
    end

    saveFilename = [saveFolder filesep sprintf('%s %s %s %s stats.xlsx',dataSet,protocol,experiment,dataType)];
    stats = myMultipleIndependentGroupStats(allData,cell_types,comments,'linear',saveFilename,savenum);

    if nCellTypes == 1
        stats.kw = [];
    end

    stats.kw.sig = sig;    
    y = getYparams(dataType,protocol,'normal');

%     if strcmp(dataType,'power')
%         y = getYparams(dataType,protocol,'log');
%     else
%         y = getYparams(dataType,protocol,'normal');
%     end
    
    groupLabels = myGroupLabels(cell_types);
    [colors,alphas] = getColorAlpha(cell_types,1);
    pTitle = sprintf('%s %s %s',dataSet,protocol,experiment);
    if strcmp(dataType,{'totalGammaPower'})
        pTitle = sprintf('%s %s %s %d-%d Hz',dataSet,protocol,experiment,gammaRange(1),gammaRange(2));
    elseif strcmp(dataType,{'totalFastPower'})
        pTitle = sprintf('%s %s %s %d-%d Hz',dataSet,protocol,experiment,fastRange(1),fastRange(2));
    end
    
    %     figBox = myBoxplot(allData,groupLabels,y,colors,alphas,stats.kw);
    figBox = myAlternativeBoxplot(allData,groupLabels,y,colors,alphas,stats.kw);
    sgtitle(pTitle,'Fontweight','bold');
    if strcmp(dataType,{'power'})
        ylim([1 5]);
    elseif any(strcmp(dataType,{'totalGammaPower','totalFastPower'}))
        if strcmpi(wavelet_norm,'L1')
            ylim([4 9]); 
        elseif strcmpi(wavelet_norm,'L2')
            ylim([2 7]); 
        end
    elseif strcmp(dataType,'frequency')
        ylim([50 150]);
    end
    
    figViolin = plotViolin(allData,groupLabels,y,colors,alphas,stats.kw);
    sgtitle(pTitle,'Fontweight','bold');
    if strcmp(dataType,{'power'})
        ylim([1 5]);
    elseif any(strcmp(dataType,{'totalGammaPower','totalFastPower'}))
        if strcmpi(wavelet_norm,'L1')
            ylim([4 9]); 
        elseif strcmpi(wavelet_norm,'L2')
            ylim([2 7]); 
        end
    elseif strcmp(dataType,'frequency')
        ylim([50 150]);
    end

    if strcmp(dataType,'power') || strcmp(dataType,'frequency')
        figBar = myBarChart(allData,groupLabels,y,colors,alphas,stats.kw);
        sgtitle(pTitle,'Fontweight','bold');
    end
    
    if savenum
        saveFilename = [saveFolder filesep sprintf('%s %s %s %s',dataSet,protocol,experiment,dataType)];
        print(figBox,'-vector','-dsvg',[saveFilename ' Boxplot.svg']); % save boxplot as svg file
        if strcmp(dataType,'power') || strcmp(dataType,'frequency')
            print(figBar,'-vector','-dsvg',[saveFilename ' BarChart.svg']); % save barchart as svg file
        end
        print(figViolin,'-vector','-dsvg',[saveFilename ' Violin.svg']); % save violin plots as svg file
    end
    
    close all;

end


end


end



end


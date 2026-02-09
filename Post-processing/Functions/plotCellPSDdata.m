function plotCellPSDdata(dataSets, comment_group1, savenum)
% independent group plots and comparisons

% dataTypes = {'PSDGammaPower','PSDFastPower','PSDfrequency'};
dataTypes = {'PSDGammaPower','PSDFastPower'};

gammaRange = [50 100];
fastRange = [100 200];

params.comments = comment_group1;
comments = cell2mat(params.comments);

% find ID parameters for analysis
params.locations = 'all';
params.cell_nums = 'all';
params.protocols = {'theta','theta_2chan'};
protocol = params.protocols{1};

% nProtocols = numel(params.protocols);

if strcmp(dataSets,'all'); dataSets = {'Camk2','Thy1','PV Transgenic','PV Viral'}; end

nDataTypes = length(dataTypes);

for iSet = 1:numel(dataSets)
dataSet = dataSets{iSet};

[info,~,data_path] = getInfo(dataSet);

switch dataSet
    case 'Thy1'
         saveFolder = 'C:\Users\brndn\Downloads\Thy1-ChR2\Raw Data\mEC\results\Summary';
%          experiments = {'excitation','inhibition'};
         experiments = {'inhibition'};
    case 'PV Transgenic'
         saveFolder = 'C:\Users\brndn\Downloads\PV-ChR2 Transgenic\Summary';
         experiments = {'inhibition'};
    case 'PV Viral'
         saveFolder = 'C:\Users\brndn\Downloads\PV-ChR2\Summary';
         experiments = {'inhibition'};
    case 'Camk2'
         saveFolder = 'C:\Users\brndn\Downloads\CaMK2-ChR2\Summary';
         experiments = {'excitation','inhibition'};
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
%         cell_types = {'fast spiking'};
    case 'excitation'
        cell_types = {'stellate','pyramidal','fast spiking'};
end
nCellTypes = length(cell_types);

allGammaPower = cell(nCellTypes,1);
allFastPower = cell(nCellTypes,1);
allFreq = cell(nCellTypes,1);

for iCell = 1:nCellTypes
    cellType = cell_types{iCell};
    params.cell_types = cell_types(iCell);

%     comments = cell2mat(params.comments);

    IDs = getIDs(info,params);

    IDs = removeIDs(IDs,info,dataSet); % skip bad recordings - get params for cells to skip
    
    if (isempty(IDs)); fprintf('No %s %s Files Found.',params.cell_types,comments); continue; end
    
    nIDs = numel(IDs);

    gammaPowerArray = cell(nIDs,1);
    fastPowerArray = cell(nIDs,1);
    freqArray = cell(nIDs,1);

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

        % compute and plot PSD from each theta cycle without artifacts 
%         [file.cell,~] = plotCyclePSD(file.cell,file.Fs,file.cycle_start_index_noArtifacts,file.cycle_length,'pwelch',file.cell.data_units,0);
%         [file.cell.psd,~] = plotCyclePSD(file.cell.raw_data,file.Fs,file.cycle_start_index_noArtifacts,file.cycle_length,'pwelch',file.cell.data_units,1);

%         freq = file.cell.psd.max_psd_gamma_freq;
%         power = file.cell.psd.sum_gamma_power;

        % compute and plot PSD from all stim data
%         [file.cell,~] =
%         plotallCyclePSD(file.cell,file.Fs,file.cycle_start_index_noArtifacts,file.cycle_length,'pwelch',file.cell.data_units,0); % old
        [file.cell.psd,~] = plotallCyclePSD(file.cell.raw_data,file.Fs,file.cycle_start_index_noArtifacts,file.cycle_length,'pwelch',file.cell.data_units,0);
%         freq = file.cell.psd.max_psd_gamma_freq_all;
%         power = file.cell.psd.sum_gamma_power_all;

        powerGamma = calculateTotalPower(file.cell.psd.S_all,file.cell.psd.f_all,gammaRange); % total gamma power: 60-100 Hz 
        powerFast = calculateTotalPower(file.cell.psd.S_all,file.cell.psd.f_all,fastRange); % total fast power: 100-200 Hz 

        freq = []; % compute autocorrelation then find peak in frequency spectra

% %       % remove low power data 
%         power_threshold = 20; % pA^2
%         removeFlag = (power > power_threshold);
%         removeFlag = ~removeFlag;
%         if removeFlag; freq = []; end

        gammaPowerArray{iID} = log10(powerGamma);
        fastPowerArray{iID} = log10(powerFast);
        freqArray{iID} = freq;
        clearvars file;
    end

    allGammaPower{iCell,1} = cell2mat(gammaPowerArray); % save cell data
    allFastPower{iCell,1} = cell2mat(fastPowerArray); % save cell data
    allFreq{iCell,1} = cell2mat(freqArray); % save cell data
end

for iType = 1:nDataTypes
    dataType = dataTypes{iType};
    
    switch dataType
        case 'PSDGammaPower'
            allData = allGammaPower;
            range = gammaRange;
        case 'PSDFastPower'
            allData = allFastPower;
            range = fastRange;
        case 'PSDfrequency'
            allData = allFreq;
    end
    
    % convert data to log for power plots and stats
    if strcmp(dataType,{'PSDGammaPower','PSDFastPower'})
        allData = cellfun(@log10,allData,'UniformOutput',false); % returns cell array
    %     saveFilename = [saveFolder filesep sprintf('%s %s %s %s LOG stats.xlsx',dataSet,protocol,experiment,dataType)];
    %     LOGstats = myMultipleIndependentGroupStats(allData,cell_types,comments,saveFilename,savenum);
    end
    
    saveFilename = [saveFolder filesep sprintf('%s %s %s %s stats.xlsx',dataSet,protocol,experiment,dataType)];
    stats = myMultipleIndependentGroupStats(allData,cell_types,comments,'linear',saveFilename,savenum);

    stats.kw.sig = 'symbol';
    
    groupLabels = myGroupLabels(cell_types);
    [colors,alphas] = getColorAlpha(cell_types,1);
    y = getYparams(dataType,protocol,'normal');
    pTitle = sprintf('%s %s %s %d-%d Hz',dataSet,protocol,experiment,range);
    
    figBox = myBoxplot(allData,groupLabels,y,colors,alphas,stats.kw);
    sgtitle(pTitle,'Fontweight','bold');
    
    figBar = myBarChart(allData,groupLabels,y,colors,alphas,stats.kw);
    sgtitle(pTitle,'Fontweight','bold');
    
    figViolin = plotViolin(allData,groupLabels,y,colors,alphas,stats.kw);
    sgtitle(pTitle,'Fontweight','bold');
    
    if savenum
        saveFilename = [saveFolder filesep sprintf('%s %s %s %s',dataSet,protocol,experiment,dataType)];
        print(figBox,'-vector','-dsvg',[saveFilename ' Boxplot.svg']); % save boxplot as svg file
        print(figBar,'-vector','-dsvg',[saveFilename ' BarChart.svg']); % save barchart as svg file
        print(figViolin,'-vector','-dsvg',[saveFilename ' Violin.svg']); % save violin plots as svg file
    end

end


end


end


    function totalPower = calculateTotalPower(S,freq,range)
        mask = (freq >= range(1) & freq <= range(2));
        % integrate PSD for total power
        df = freq(2) - freq(1);
        totalPower = trapz(S(mask))*df;
    end



end


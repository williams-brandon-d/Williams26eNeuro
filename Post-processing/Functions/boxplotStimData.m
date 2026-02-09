function boxplotStimData(dataSets, comment_group1, savenum)
% independent group plots and comparisons

% dataTypes = {'stim','power','frequency','phase'};
dataTypes = {'stim'};

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
         experiments = {'excitation','inhibition','currentclamp'};
%          experiments = {'inhibition'};
    case 'PV Transgenic'
         saveFolder = 'C:\Users\brndn\Downloads\PV-ChR2 Transgenic\Summary';
         experiments = {'inhibition','currentclamp'};
    case 'PV Viral'
         saveFolder = 'C:\Users\brndn\Downloads\PV-ChR2\Summary';
         experiments = {'inhibition','currentclamp'};
    case 'Camk2'
         saveFolder = 'C:\Users\brndn\Downloads\CaMK2-ChR2\Summary';
         experiments = {'excitation','inhibition','currentclamp'};
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

allStim = cell(nCellTypes,1);

for iCell = 1:nCellTypes
    cellType = cell_types{iCell};
    params.cell_types = cell_types(iCell);

    IDs = getIDs(info,params);

    IDs = removeIDs(IDs,info,dataSet); % skip bad recordings - get params for cells to skip
    
    if (isempty(IDs)); fprintf('No %s %s Files Found.',params.cell_types,comments); continue; end
    
    nIDs = numel(IDs);

    stimArray = cell(nIDs,1);

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
        stim = file.led_input;
        % convert stim led input to light intensity

        stimArray{iID} = stim;

        clearvars file;
    end

    allStim{iCell,1} = cell2mat(stimArray); % save cell data

end


for iType = 1:nDataTypes
    dataType = dataTypes{iType};
    
    switch dataType
        case 'stim'
            allData = allStim;
    end

    saveFilename = [saveFolder filesep sprintf('%s %s %s %s stats.xlsx',dataSet,protocol,experiment,dataType)];
    stats = myMultipleIndependentGroupStats(allData,cell_types,comments,'linear',saveFilename,savenum);

    if nCellTypes == 1
        stats.kw = [];
    end

    stats.kw.sig = sig;
    
    y = getYparams(dataType,protocol,'normal');
    
    groupLabels = myGroupLabels(cell_types);
    [colors,alphas] = getColorAlpha(cell_types,1);
    pTitle = sprintf('%s %s %s',dataSet,protocol,experiment);
    
%     figBox = myBoxplot(allData,groupLabels,y,colors,alphas,stats.kw);
    figBox = myAlternativeBoxplot(allData,groupLabels,y,colors,alphas,stats.kw);
    sgtitle(pTitle,'Fontweight','bold');
    
    figViolin = plotViolin(allData,groupLabels,y,colors,alphas,stats.kw);
    sgtitle(pTitle,'Fontweight','bold');
    
    
    if savenum
        saveFilename = [saveFolder filesep sprintf('%s %s %s %s',dataSet,protocol,experiment,dataType)];
        print(figBox,'-vector','-dsvg',[saveFilename ' Boxplot.svg']); % save boxplot as svg file
        print(figViolin,'-vector','-dsvg',[saveFilename ' Violin.svg']); % save violin plots as svg file
    end
    
    close all;

end


end


end



end


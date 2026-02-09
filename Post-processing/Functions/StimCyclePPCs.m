function StimCyclePPCs(dataSets, savenum)

titlefontsize = 15;
cell_types = {'stellate','pyramidal','fast spiking'};

% find ID parameters for analysis
params.locations = 'all';
params.experiments = {'currentclamp'};
params.cell_nums = 'all';
params.comments = {'','DNQX before','GBZ before' 'DNQX_GABAzine_before', 'Gabazine before',}; % 'GBZ before' 'DNQX_GABAzine_before', 'Gabazine before', 'DNQX before'
params.protocols = {'theta','theta_2chan'};

if strcmp(dataSets,'all'); dataSets = {'Camk2','Thy1','PV Transgenic','PV Viral'}; end

% load data for figure

nCellTypes = numel(cell_types);

for iSet = 1:numel(dataSets)
dataSet = dataSets{iSet};

[info,~,data_path] = getInfo(dataSet);

switch dataSet
    case 'Thy1'
         saveFolder = 'C:\Users\brndn\Downloads\Thy1-ChR2\Raw Data\mEC\results\Summary';
    case 'PV Transgenic'
         saveFolder = 'C:\Users\brndn\Downloads\PV-ChR2 Transgenic\Summary';
    case 'PV Viral'
         saveFolder = 'C:\Users\brndn\Downloads\PV-ChR2\Summary';
    case 'Camk2'
         saveFolder = 'C:\Users\brndn\Downloads\CaMK2-ChR2\Summary';
end

if ~exist(saveFolder, 'dir')
   mkdir(saveFolder)
end

% allData = cell(nCellTypes,1);
% allColor = cell(nCellTypes,1);
% allN = cell(nCellTypes,1);
% allnCycles = cell(nCellTypes,1);

for iCell = 1:nCellTypes
    cellType = cell_types{iCell};
    params.cell_types = cell_types(iCell);
    IDs = getIDs(info,params);
    
    if (isempty(IDs)); fprintf('No %s Files Found.',params.cell_types); continue; end
    
    nIDs = numel(IDs);

    dataCell = cell(nIDs,1);
    
    for iID = 1:nIDs
    
        ID = IDs{iID};
        ID_index = find_index(info,'ID',ID);
        filename = sprintf('%s.abf',ID);
        fprintf('Analyzing %s,File %d/%d\n',filename,iID,nIDs)
        p = info(ID_index);
    
        if isempty(p.comments)
            comment = 'No Comment';
        else
            comment = p.comments;
        end
    
        commentFolder = sprintf('%sresults\\%s\\%s\\%s\\%s\\%s\\%s\\%s',data_path,p.location,p.cell_type,p.cell_num,p.experiment,p.protocol,comment);
    
        dataFolder = getIDFolder(commentFolder,ID);
    
        dataFilename = [dataFolder filesep 'data.mat'];
    
        load(dataFilename,'file');

        % gather relevant data 
        dataCell{iID} = file.theta_spike_phases; % nCycles x 1 cell array
    
        clearvars file;
    end
    
    % plot phase locking for each cell type
    figCell = plotCyclePPCs(dataCell,cellType); % plot figure for each cell type
    titleString = {sprintf('%s %s',dataSet,cellType);'Theta Stim Phase Locking'};
    title(titleString,'FontSize',titlefontsize,'Fontweight','bold');
    
    if savenum
        saveFilename = [saveFolder filesep sprintf('%s %s Theta Stim Phase Locking.svg',dataSet,cellType)];
        print(figCell,'-vector','-dsvg',saveFilename);
    end
    
    close all;
% 
%     allData{iCell} = [dataCell{:}]; % save all cell data
%     allColor{iCell} = color;
%     allN{iCell} = nIDs;
%     allnCycles{iCell} = nCycles;

end

% figAll = plotAvgHist(allData,allN,allnCycles,allColor,cell_types);
% pTitle = sprintf('%s Theta Stim Phase Histogram Average',dataSet);
% sgtitle(pTitle,'Fontweight','bold');
% 
% if savenum
%     saveFilename = [saveFolder filesep sprintf('%s Theta Stim Phase Histogram Average.svg',dataSet)];
%     print(figAll,'-vector','-dsvg',saveFilename);
% end

end



end
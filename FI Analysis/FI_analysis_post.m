%% compare features across excitatory cells

% Outliers
% stellate id 21819031 has long latency
% stellate id 20901000 and 20910000 have long time constant

% compare mean ISI ratio from all current levels

% combine thy1 and pv datasets?

clear variables; close all; clc;
cd('C:\Users\brndn\OneDrive\Desktop\White Lab\Matlab Files\');
addpath(genpath('./')); % add all folders and subfolders in cd to path

savenum = 1; % 0 = don't save, 1 = save info

dataSets = {'Thy1'};

params.protocols = {'FI'};
protocol = params.protocols{1};
params.experiments = {'currentclamp'};

params.locations = 'all';
params.cell_nums = 'all';
 
cell_types = {'stellate','pyramidal'};
params.comments = {'','DNQX before','Gabazine before'}; % combine comment data

% features = {'ISI_ratio','sag_hyper','sag_depol_thresh','dap','latency','time_constant','rmp'};
features = {'rmp'};

sig = 'exact';

comments = cell2mat(params.comments);

nFeatures = length(features);
nCellTypes = length(cell_types);

tic;
for iSet = 1:numel(dataSets)
dataSet = dataSets{iSet};

[info,info_fullname,data_path] = getInfo(dataSet);

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

allISI = cell(nCellTypes,1);
allsagHyper = allISI;
allsagDepol = allISI;
alldap = allISI;
alllatency = allISI;
alltau = allISI;
allrmp = allISI;

for iCell = 1:nCellTypes
    cellType = cell_types{iCell};
    params.cell_types = cell_types(iCell);

    IDs = getIDs(info,params);
    
    if (isempty(IDs)); fprintf('No %s %s Files Found.',params.cell_types,comments); continue; end
    
    nIDs = numel(IDs);

    ISIArray = cell(nIDs,1);
    sagHyperArray = ISIArray;
    sagDepolArray = ISIArray;
    dapArray = ISIArray;
    latencyArray = ISIArray;
    tauArray = ISIArray;
    rmpArray = ISIArray;

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
        ISIArray{iID} = file.ISI_ratio;
        sagHyperArray{iID} = file.sag_hyper;
        sagDepolArray{iID} = file.sag_depol_thresh;
        dapArray{iID} = file.dap;
        latencyArray{iID} = file.latency;
        tauArray{iID} = file.time_constant;
        rmpArray{iID} = file.RMP;

        clearvars file;
    end

    allISI{iCell,1} = cell2mat(ISIArray); % save cell data
    allsagHyper{iCell,1} = cell2mat(sagHyperArray); % save cell data
    allsagDepol{iCell,1} = cell2mat(sagDepolArray); % save cell data
    alldap{iCell,1} = cell2mat(dapArray); % save cell data
    alllatency{iCell,1} = cell2mat(latencyArray); % save cell data
    alltau{iCell,1} = cell2mat(tauArray); % save cell data
    allrmp{iCell,1} = cell2mat(rmpArray); % save cell data
end

for iFeature = 1:nFeatures
    dataType = features{iFeature};
    
    switch dataType
        case 'ISI_ratio'
            allData = allISI;
        case 'sag_hyper'
            allData = allsagHyper;
        case 'sag_depol_thresh'
            allData = allsagDepol;
        case 'dap'
            allData = alldap;
        case 'latency'
            allData = alllatency;
        case 'time_constant'
            allData = alltau;
        case 'rmp'
            allData = allrmp;
    end
    
    saveFilename = [saveFolder filesep sprintf('%s %s %s stats.xlsx',dataSet,protocol,dataType)];
    stats = myMultipleIndependentGroupStats(allData',cell_types,comments,'linear',saveFilename,savenum);
    
    groupLabels = myGroupLabels(cell_types);
    [colors,alphas] = getColorAlpha(cell_types,1);
    y = getYparams(dataType,protocol,'normal');
    pTitle = sprintf('%s',dataSet);

    stats.rs.sig = sig;
    
    figBox = myBoxplot(allData,groupLabels,y,colors,alphas,stats.rs);
    sgtitle(pTitle,'Fontweight','bold');
    
    figBar = myBarChart(allData,groupLabels,y,colors,alphas,stats.rs);
    sgtitle(pTitle,'Fontweight','bold');
    
    figViolin = plotViolin(allData,groupLabels,y,colors,alphas,stats.rs);
    sgtitle(pTitle,'Fontweight','bold');
    
    if savenum
        saveFilename = [saveFolder filesep sprintf('%s %s %s',dataSet,protocol,dataType)];
        print(figBox,'-vector','-dsvg',[saveFilename ' Boxplot.svg']); % save boxplot as svg file
        print(figBar,'-vector','-dsvg',[saveFilename ' BarChart.svg']); % save barchart as svg file
        print(figViolin,'-vector','-dsvg',[saveFilename ' Violin.svg']); % save violin plots as svg file
    end

end


end

function plotInterspikeFiringRates(dataSets, savenum)
% 

cell_types = {'stellate','pyramidal','fast spiking'};

% find ID parameters for analysis
params.locations = 'all';
params.experiments = {'currentclamp'};
params.cell_nums = 'all';
params.comments = {'','DNQX before','DNQX before 10 uM','GBZ before' 'DNQX_GABAzine_before', 'Gabazine before'}; % 'GBZ before' 'DNQX_GABAzine_before', 'Gabazine before', 'DNQX before'
params.protocols = {'theta','theta_2chan'};

if strcmp(dataSets,'all'); dataSets = {'Camk2','Thy1','PV Transgenic','PV Viral'}; end

comments = cell2mat(params.comments);

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

allData = cell(nCellTypes,1);
allColor = cell(nCellTypes,1);

for iCell = 1:nCellTypes
    cellType = cell_types{iCell};
    params.cell_types = cell_types(iCell);
    IDs = getIDs(info,params);
    
    if (isempty(IDs)); fprintf('No %s Files Found.',params.cell_types); continue; end
    
    nIDs = numel(IDs);
    
    switch cellType
        case 'stellate'
            color = [1 0 0];
        case 'pyramidal'
            color = [1 0.4 0];
        case 'fast spiking'
            color = [0 0 1];
    end

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

        % get cell FRs
        dataCell{iID} = cell2mat(file.FRs(2:end)); % skip first stim period
        
        clearvars file;
    end
    
    allData{iCell} = vertcat(dataCell{:}); % save all cell data
    allColor{iCell} = color;

end

% combine stellate and pyramidal cells?

% data table and stats
saveFilename = [saveFolder filesep sprintf('%s Interspike Firing Rate stats.xlsx',dataSet)];
stats = myMultipleIndependentGroupStats(allData,cell_types,comments,'linear',saveFilename,savenum);

% alternative boxplots or violin?

fig = plotViolin(allData,allColor,cell_types);
pTitle = sprintf('%s',dataSet);
sgtitle(pTitle,'Fontweight','bold');

if savenum
    saveFilename = [saveFolder filesep sprintf('%s Interspike Firing Rate Average.svg',dataSet)];
    print(fig,'-vector','-dsvg',saveFilename);
end

end

    function fig = plotViolin(allData,allColor,cell_types)
        % add nice colors
%         c =  [0.45, 0.80, 0.69;...
%               0.98, 0.40, 0.35;...
%               0.55, 0.60, 0.79;...
%               0.90, 0.70, 0.30]; 

        % check if any columns are all NaN - remove group but save correct labels
        nGroups = numel(allData);

        % pad with NaN
        maxNumEl = max(cellfun(@numel,allData(:)));
        data_all_pad = cellfun(@(x){padarray(x(:),[maxNumEl-numel(x),0],NaN,'post')}, allData(:)); % Pad each vector with NaN values to equate lengths
        dataArray = cell2mat(data_all_pad'); 

        linData = dataArray(:);

        nan_mask = isnan(linData);
        linData(nan_mask) = []; % remove NaNs from linear array

        groups = ones(numel(dataArray),1);
        for i = 2:nGroups
            idx_start = ((i-1)*size(dataArray,1)) + 1;
            idx_stop = idx_start + size(dataArray,1) - 1;
            groups(idx_start:idx_stop) = i*ones(size(dataArray,1),1);
        end
        groups(nan_mask) = []; % remove NaNs from linear array

        newGroups = unique(groups);
        newNGroups = numel(newGroups);

        colors = cell2mat(allColor(newGroups));

%         labels = cell(newNGroups,1);
        row1 = cell(1,newNGroups);
        row2 = row1;
        for iGroup = 1:newNGroups
            group = newGroups(iGroup);
            data_group = linData(groups == group);
            N = numel(data_group);
            row1{iGroup} = sprintf('%s',cell_types{iGroup});
            row2{iGroup} = sprintf('(n=%d)',N);
        end
        labelArray = [row1;row2];
        tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));

        tickfontsize = 15;

        % plot violin
        fig = figure;
        if ~isempty(linData)
            if nargin < 3 
                daviolinplot(linData,'groups',groups); % default colors
            else
                daviolinplot(linData,'groups',groups,'color',colors); % my colors
            end
        end
        ylim([0 inf])
        ylabel({'Avg Theta Firing Rate'});
        ax = gca;
        ax.TickLabelInterpreter = 'tex';
        ax.XTickLabel = tickLabels;
        ax.YAxis.FontSize = tickfontsize;
        ax.YAxis.FontWeight = 'bold';
        ax.XAxis.FontSize = tickfontsize;
        ax.XAxis.FontWeight = 'bold';
    end



end
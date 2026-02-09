function Thy1_DorsalvsVentral()
clear variables; close all; clc;
cd('C:\Users\brndn\OneDrive\Desktop\White Lab\Matlab Files\');

% there are more recordings for pulse inhibition vs theta
% only 1 FS for dorsal and ventral for pulse inhibition - 0 for theta

savenum = 1; % 0 = don't save, 1 = save info
dataSet = 'Thy1';
saveFolder = 'C:\Users\brndn\Downloads\Thy1-ChR2\Raw Data\mEC\results\Summary';

% experiments = {'inhibition','excitation','currentclamp'};
dataTypes = {'power','frequency','rate'};
cell_types = {'stellate','pyramidal'};
locations = {'dorsal','ventral'};

% find ID parameters for analysis
params.cell_nums = 'all';
params.comments = {'','DNQX before','DNQX before 10 uM','GBZ before' 'DNQX_GABAzine_before', 'Gabazine before'}; % normal conditions
params.protocols = {'theta'};
protocol = params.protocols{1};

sig = 'symbol'; % 'symbol' or 'exact' | p-value significance indicator


comments = cell2mat(params.comments);

% nExperiments = numel(experiments);
nDataTypes = numel(dataTypes);
nCellTypes = numel(cell_types);
nLocations = numel(locations);

[info,~,data_path] = getInfo(dataSet);

if ~exist(saveFolder, 'dir')
   mkdir(saveFolder)
end


% for each dataType
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

allData = cell(nCellTypes,nLocations);
allColor = cell(nCellTypes,nLocations);
allAlpha = cell(nCellTypes,nLocations);

for iCell = 1:nCellTypes
    cellType = cell_types{iCell};
    params.cell_types = cell_types(iCell);

    switch cellType
        case 'stellate'
            color = [1 0 0];
        case 'pyramidal'
            color = [1 0.4 0];
        case 'fast spiking'
            color = [0 0 1];
    end


    for iLoc = 1:nLocations
        location = locations{iLoc};
        params.locations = {location};

        switch location
            case 'dorsal'
                alpha = 1;
            case 'ventral'
                alpha = 0.4;
        end

        IDs = getIDs(info,params);

        % check for any bad recordings?
        IDs = removeIDs(IDs,info,dataSet); % skip bad recordings - get params for cells to skip

        if (isempty(IDs)); fprintf('No %s Files Found.',cellType); continue; end
        
        nIDs = numel(IDs);
    
        data = zeros(nIDs,1);
        
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
            switch dataType 
                case 'power'
                    if ~isempty(file.cell.CWTmaxValues)
                        data(iID) = log10(file.cell.CWTmaxValues(3)); % log lfp peak scalogram power
                    else
                        data(iID) = NaN;
                    end
                case 'frequency'
                    if ~isempty(file.cell.CWTmaxValues)
                        % if data is frequency or phase - remove low power data 
                        removeFlag = removeLowPowerData(file,dataSet);
                        if ~removeFlag
                            data(iID) = file.cell.CWTmaxValues(2); % lfp peak scalogram freq
                        else
                            data(iID) = []; 
                        end
%                         data(iID) = file.cell.CWTmaxValues(2); % lfp peak scalogram freq
                    else
                        data(iID) = NaN;
                    end
                case 'rate'
                    % calculate average firing rate
    %                 file.FRs = FRs;
                    nSpikes = numel([file.theta_spike_phases{:}]); % nCycles x cycle phase array
                    nCycles = numel(file.cycles);
                    if ~isempty(nSpikes)
                        data(iID) = nSpikes / nCycles; % average spikes per theta cycle 
                    else
                        data(iID) = 0; % change to NaN to ignore non-firing cells
                    end
            end

            clearvars file;

        end

        allData{iCell,iLoc} = data; % save all cell data
        allColor{iCell,iLoc} = color;
        allAlpha{iCell,iLoc} = alpha;

    end

end

% data table and stats
saveFilename = [saveFolder filesep sprintf('%s Dorsal vs Ventral %s stats.xlsx',dataSet,dataType)];
stats = myMultipleIndependentGroupStats(allData,cell_types,comments,'linear',saveFilename,savenum);

% combine E cells stats?
combinedData = cell(1,nLocations);
for iLoc = 1:nLocations
    combinedData{iLoc} = vertcat(allData{:,iLoc});
end

% save combined E cell stats
saveFilenameCombined = [saveFolder filesep sprintf('%s Dorsal vs Ventral %s combined E cells stats.xlsx',dataSet,dataType)];
combinedStats = myMultipleIndependentGroupStats(combinedData,{'E cells'},locations,'linear',saveFilenameCombined,savenum);

stats.kw.sig = sig;
combinedStats.kw.sig = sig;

% boxplots with all cell types and dorsal/ventral

switch dataType
    case {'power','frequency'}
        y = getYparams(dataType,protocol,'normal');
    case 'rate'
    y.min = 0; y.max = 5; y.dy = 1;
    y.ticks = y.min:y.dy:y.max;
    y.tickLabels = compose('%d',y.ticks);
    y.labelstring = 'Firing Rate (Spikes / Theta Cycle)';   
    y.scale = 'linear';
end

% groupLabels = myGroupLabels(cell_types);
% [colors,alphas] = getColorAlpha(cell_types,1);
pTitle = sprintf('%s %s %s',dataSet,protocol,experiment);

% fig = plotViolin(allData,allColor,cell_types);
fig = myAlternativeBoxplot(allData,cell_types,y,allColor,allAlpha,stats);
if strcmp(dataType,'frequency')
    ylim([50 150]);
end
sgtitle(pTitle,'Fontweight','bold');

figCombined = myAlternativeBoxplot(combinedData,{'E cells'},y,allColor,allAlpha,combinedStats);
if strcmp(dataType,'frequency')
    ylim([50 150]);
end
sgtitle(pTitle,'Fontweight','bold');

if savenum
    saveFilename = [saveFolder filesep sprintf('%s Dorsal vs Ventral %s .svg',dataSet,dataType)];
    print(fig,'-vector','-dsvg',saveFilename);
    saveFilenameCombined = [saveFolder filesep sprintf('%s Dorsal vs Ventral %s combined .svg',dataSet,dataType)];
    print(figCombined,'-vector','-dsvg',saveFilenameCombined);
end


end

    function fig = myAlternativeBoxplot(data_all,xTickLabels,y,colors,alphas,stats)
     % Input: data_all - nCellTypes x nGroups

    colors = cell2mat(colors);
    
    tickfontsize = 15;
    scattersize = 70;
    
    [nCells,nGroups] = size(data_all);
    
    xPos = 1:nCells;
    
    data_trans = data_all'; % for easier grouping transpose so groups are in order
    maxNumEl = max(cellfun(@numel,data_trans(:)));
    data_all_pad = cellfun(@(x){padarray(x(:),[maxNumEl-numel(x),0],NaN,'post')}, data_trans); % Pad each vector with NaN values to equate lengths
    data_all_mat = cell2mat(data_all_pad); 
    
    data_groups = zeros(size(data_all_mat,1),1);
    count = 0;
    for iii = 1:nGroups
        count = count + 1;
        ix1 = (count-1).*(maxNumEl) + 1;
        ix2 = ix1 + maxNumEl - 1;
        data_groups(ix1:ix2) = count.*ones(maxNumEl,1);
    end

    % data_all_mat needs to be 3 columns containing both dorsal ventral data points
    % data_groups must have group indices corresponding to data_all_mat columns

    RGB = crameri('roma',nGroups);

    groups = {'dorsal','ventral'};
    
    
    fig = figure;
    
    if ~isempty(data_all_mat)
    %     plot_handle = boxplot(data_all_mat,data_groups,'Colors','k','Symbol','k.','OutlierSize',10,'Positions',xPos);
    
    %     h = daboxplot(data_all_mat,'groups',data_groups,'xtlabels', xTickLabels,...
    %         'colors',colors,'boxalpha',0.7,'whiskers',0,'scatter',1,'outsymbol','k*',...
    %         'outliers',1,'scattersize',25,'flipcolors',0,'boxspacing',1.2,'scattercolors',{'k','w'}); 
    
        h = daboxplot(data_all_mat,'groups',data_groups,'xtlabels', xTickLabels,...
            'colors',RGB,'fill',1,'boxalpha',0.4,'whiskers',0,'scatter',1,'outsymbol','k*',...
            'outliers',1,'scattersize',scattersize,'flipcolors',0,'boxspacing',1.2,'scattercolors',{'k','w'},'legend',groups); 
    
%         for g = 1:numel(xTickLabels)
%            set(h.sc(:,g),'MarkerFaceColor',colors(g,:));
%            set(h.ot(:,g),'CData',colors(g,:));
%         end
    
    
    end
    
    % box off

    legend('boxoff');

    % xticks(xPos)
    % xticklabels(xTickLabels)

%     han = legend(xTickLabels);
    
    ax = gca(); 
    ax.TickLabelInterpreter = 'tex';  % needed for some plots like boxplot.
    
    ax.XAxis.FontSize = tickfontsize;
    ax.XAxis.FontWeight = 'bold';
    % xtickangle(ax,45)
    
    ylabel(y.labelstring,'FontSize',25,'FontWeight','bold')
    ax.YLim = [y.min y.max];
    ax.YAxis.Scale = y.scale; 
    ax.YTick = y.ticks;
    ax.YTickLabel = y.tickLabels;
    ax.YAxis.FontSize = tickfontsize;
    ax.YAxis.FontWeight = 'bold';
    
    if isfield(stats,'results') 
        if isfield(stats,'tbl')
            pvalue = stats.tbl{1,end};
            pvalue = pvalue{1};
            if pvalue > 0.05 % check pvalue for anova or kw test
                return;
            end
        end
            
        dataMax = max(data_all_mat(:));
        plotSignificance(stats,xPos,dataMax,y.scale)
    end

end



end
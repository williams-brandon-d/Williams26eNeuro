function plotThetaAmplitudeValues(dataSets, comment_group1, savenum)
% independent group plots and comparisons

dataTypes = {'theta amplitude'};

params.comments = comment_group1;
comments = cell2mat(params.comments);

% find ID parameters for analysis
params.locations = 'all';
params.cell_nums = 'all';
params.protocols = {'theta'};
% protocol = params.protocols{1};

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
    case 'excitation'
        cell_types = {'stellate','pyramidal','fast spiking'};
end
nCellTypes = length(cell_types);

for iCell = 1:nCellTypes
    cellType = cell_types{iCell};
    params.cell_types = cell_types(iCell);

%     comments = cell2mat(params.comments);

    IDs = getIDs(info,params);

    IDs = removeIDs(IDs,info); % skip bad recordings - get params for cells to skip
    
    if (isempty(IDs)); fprintf('No %s %s Files Found.',params.cell_types,comments); continue; end
    
    nIDs = numel(IDs);

    phaseArray = cell(nIDs,1);
    freqArray = cell(nIDs,1);
    powerArray = cell(nIDs,1);

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
        if ~isfield(file.cell,'CWTcycleValues')
            [~,~,~,file.cell.CWTcycleValues] = thetaCWT(file,'cell');
            save(dataFilename,'file'); % save data
        end

        % if cycle has artifact - assign NaN value
        nCycles = numel(file.cycles);
        cycleValues = NaN(3,nCycles);
        iNoArtifact = 0;
        for iCycle = 1:nCycles
            if ~file.cell.spikes(iCycle)
                iNoArtifact = iNoArtifact + 1;
                cycleValues(:,iCycle) = file.cell.CWTcycleValues(iNoArtifact,:);
            end
        end

        phaseArray{iID} = cycleValues(1,:);
        freqArray{iID} = cycleValues(2,:);
        powerArray{iID} = cycleValues(3,:);

        cycles = file.cycles; % all data should be processed for the same cycles
        clearvars file;
    end

    [color,~] = getCellcolor(cellType,dataSet,comment);
%     exampleHist = getExamplePhaseHist(exampleID,info,data_path,edges);

    for iType = 1:nDataTypes
        dataType = dataTypes{iType};

        switch dataType
            case 'power'
                data = log10(cell2mat(powerArray)); % log transform power data
            case 'frequency'
                data = cell2mat(freqArray);
            case 'phase'
                data = cell2mat(phaseArray);
        end
        
        figCell = plotCWTcycleValues(data,cycles,color,dataType); % plot figure for each cell type
    
        pTitle = sprintf('%s %s',dataSet,cellType);
        sgtitle(pTitle,'Fontweight','bold');
        
        if savenum
            saveFilename = [saveFolder filesep sprintf('%s %s Theta CWT Peak %s over time.svg',dataSet,cellType,dataType)];
            print(figCell,'-vector','-dsvg',saveFilename);
        end

    end
    
%     close all;

end

end

end

end

function fig = plotCWTcycleValues(data,xPos,color,dataType)
% data = nIDs x nCycles matrix

% calculate mean and SEM for each cycle

tickfontsize = 15;
    
fig = figure;
hold on;

if ~isempty(data)
    % plot SEM
    meanData = mean(data,1,'omitnan');
    SEM = std(data,1,'omitnan') ./ sqrt(size(data,1));
    meanData = reshape(meanData,1,[]); % column vector
    SEM = reshape(SEM,1,[]); % column vector
    xSEM = [xPos fliplr(xPos)] ;         
    ySEM = [meanData+SEM fliplr(meanData-SEM)];

    han1 = fill(xSEM,ySEM,color);
    han1.FaceColor = color;    
    han1.FaceAlpha = 0.4;      
    han1.EdgeColor = 'none'; 
    drawnow;

    % plot mean
    plot(xPos,meanData,'Color',color);

%     % plot example
%     if ~isempty(exampleData)
%         plot(xPos,exampleData,'Color',color,'LineStyle','--');
%     end

end

hold off

xlim([-inf inf])
xlabel('Theta Cycle #')

% ylim([-inf inf])
switch dataType
    case 'power'
        ylabel('Log Peak Gamma Power pA^{2}');
        yLim = [2 3.5]; dy = 0.5;
        ylim(yLim);
        yTicks = yLim(1):dy:yLim(2);
        yticks(yTicks);
        yticklabels(string(yTicks));
    case 'frequency'
        ylabel('Frequency (Hz)');
        yLim = [50 150]; dy = 50;
        ylim(yLim);
        yTicks = yLim(1):dy:yLim(2);
        yticks(yTicks);
        yticklabels(string(yTicks));
    case 'phase'
        ylabel('Theta Phase (rad)');
        yLim = [-pi/2 pi/2]; dy = pi/2;
        ylim(yLim);
        yTicks = yLim(1):dy:yLim(2);
        yticks(yTicks);
        yticklabels({'-π/2','0','π/2'});
end

ax = gca;
ax.YAxis.FontSize = tickfontsize;
ax.YAxis.FontWeight = 'bold';
ax.XAxis.FontSize = tickfontsize;
ax.XAxis.FontWeight = 'bold';


end


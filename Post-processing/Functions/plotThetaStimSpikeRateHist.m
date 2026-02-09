function plotThetaStimSpikeRateHist(dataSets, savenum)
% histogram all FRs for each cell type

edges = 0:10:400; % Hz
nbins = numel(edges) - 1;

% find ID parameters for analysis
params.locations = 'all';
params.experiments = {'currentclamp'};
params.cell_nums = 'all';
params.comments = {'','DNQX before','DNQX before 10 uM','GBZ before' 'DNQX_GABAzine_before', 'Gabazine before',}; % 'GBZ before' 'DNQX_GABAzine_before', 'Gabazine before', 'DNQX before'
params.protocols = {'theta','theta_2chan'};

if strcmp(dataSets,'all'); dataSets = {'Camk2','Thy1','PV Transgenic','PV Viral'}; end

for iSet = 1:numel(dataSets)
dataSet = dataSets{iSet};

[info,~,data_path] = getInfo(dataSet);

switch dataSet
    case 'Thy1'
         saveFolder = 'C:\Users\brndn\Downloads\Thy1-ChR2\Raw Data\mEC\results\Summary';
         cell_types = {'stellate','pyramidal','fast spiking'};
    case 'PV Transgenic'
         saveFolder = 'C:\Users\brndn\Downloads\PV-ChR2 Transgenic\Summary';
         cell_types = {'fast spiking'};
    case 'PV Viral'
         saveFolder = 'C:\Users\brndn\Downloads\PV-ChR2\Summary';
         cell_types = {'fast spiking'};
    case 'Camk2'
         saveFolder = 'C:\Users\brndn\Downloads\CaMK2-ChR2\Summary';
         cell_types = {'stellate','pyramidal','fast spiking'};
end

if ~exist(saveFolder, 'dir')
   mkdir(saveFolder)
end

nCellTypes = numel(cell_types);

allData = cell(nCellTypes,1);
colors = cell(nCellTypes,1);

for iCell = 1:nCellTypes
    cellType = cell_types{iCell};
    params.cell_types = cell_types(iCell);
    IDs = getIDs(info,params);

    switch dataSet
        case 'Thy1'
            switch cellType 
                case 'stellate'
                    exampleID = '21317023'; % stellate 16
                case 'pyramidal'
                    exampleID = '21120002'; % pyramidal 4
                case 'fast spiking'
                    exampleID = '21217013'; % fast spiking 18
            end
        case 'PV Transgenic'
            switch cellType 
                case 'fast spiking'
                    exampleID = '23228040'; % fast spiking 8
            end
        case 'Camk2'
            switch cellType
                case 'stellate'
                    exampleID = '23o12022'; % stellate 14
                case 'pyramidal'
                    exampleID = '23o11002'; % pyramidal 9
                case 'fast spiking'
                    exampleID = '23o18036'; % fast spiking 9
            end
    end
    
    if (isempty(IDs)); fprintf('No %s Files Found.',params.cell_types); continue; end
    
    nIDs = numel(IDs);

    HistArray = zeros(nIDs,nbins);
    
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
        FRs = file.FRs;
%         cell_FRs = cell2mat(FRs); % nCycles x cycle phase array

%         removeNum = 2; % num of beginning cycles to remove from analysis
%         removeIndices = false(numel(FRs),1);
%         if removeNum
%             removeIndices(1:removeNum) = 1; % add to removeIndices
%         end
%         FRs = FRs(~removeIndices);

        nCycles = numel(FRs);

        cell_FRs = cell2mat(FRs); % nCycles x cycle phase array
%         cell_FRs = [FRs{:}]; % nCycles x cycle phase array

%         nCycles = numel(file.cycles);

        if ~isempty(cell_FRs)
            [HistArray(iID,:),~] = histcounts(cell_FRs,edges);
            HistArray(iID,:) = HistArray(iID,:)/nCycles;
        end
    
        clearvars file;
    end
    
    [color,~] = getCellcolor(cellType,dataSet,comment);
    exampleHist = getExampleRateHist(exampleID,info,data_path,edges);

    fig = plotHist(HistArray,edges,color,exampleHist,dataSet); % plot figure for each cell type
    pTitle = sprintf('%s %s',dataSet,cellType);
    sgtitle(pTitle,'Fontweight','bold');
    
    if savenum
        saveFilename = [saveFolder filesep sprintf('%s %s Firing Rate Histogram.svg',dataSet,cellType)];
        print(fig,'-vector','-dsvg',saveFilename);
    end
    
%     close all;
    
    allData{iCell} = HistArray; % save all cell data
    colors{iCell} = color;

end

% stats for firing rates


figAll = plotAvgHist(allData,edges,cell_types,colors);
pTitle = sprintf('%s',dataSet);
sgtitle(pTitle,'Fontweight','bold');

if savenum
    saveFilename = [saveFolder filesep sprintf('%s Firing Rate Histogram Average.svg',dataSet)];
    print(figAll,'-vector','-dsvg',saveFilename);
end

end




    function fig = plotHist(data,edges,Color,exampleData,dataSet)

        % force edges to be a column vector
        edges = reshape(edges,[],1);

        % plot histogram of average
        delta_bin = edges(2) - edges(1);
        bin_midpoints = edges(1:end-1) + delta_bin/2;
        bin_midpoints = reshape(bin_midpoints,1,[]); % column vector

        tickfontsize = 15;
            
        fig = figure;
        hold on;
    
%         plot(bin_midpoints,HistArray,'LineStyle','--','Color',color);
%         plot(bin_midpoints,mean(HistArray,1),'LineStyle','-','Color',color,'Linewidth',2);

        if ~isempty(data)
            % plot SEM
            meanData = mean(data,1);
            SEM = std(data,1) ./ sqrt(size(data,1));
            meanData = reshape(meanData,1,[]); % column vector
            SEM = reshape(SEM,1,[]); % column vector
            xSEM = [bin_midpoints fliplr(bin_midpoints)] ;         
            ySEM = [meanData+SEM fliplr(meanData-SEM)];

            han1 = fill(xSEM,ySEM,Color);
            han1.FaceColor = Color;    
            han1.FaceAlpha = 0.4;      
            han1.EdgeColor = 'none'; 
            drawnow;

            % plot mean
            plot(bin_midpoints,meanData,'Color',Color);

            % plot example
            if ~isempty(exampleData)
                plot(bin_midpoints,exampleData,'Color',Color,'LineStyle','--');
            end


        end
    
        hold off
        xlim([min(edges) max(edges)])
        switch dataSet
            case 'Thy1'
                ylim([0 0.6]);
            case 'PV Transgenic'
                ylim([0 2]);
            case 'Camk2'
                ylim([0 0.6])
        end
        xlabel('Interspike Firing Rate (Hz)')
        ylabel('Counts (1 / Theta Cycle)')
       
        ax = gca;
        ax.YAxis.FontSize = tickfontsize;
        ax.YAxis.FontWeight = 'bold';
        ax.XAxis.FontSize = tickfontsize;
        ax.XAxis.FontWeight = 'bold';
    end

    function fig = plotAvgHist(allData,edges,cell_types,colors)

        % plot histogram of average

        edges = reshape(edges,1,[]);

        nTypes = numel(cell_types);

        delta_bin = edges(2) - edges(1);
        bin_midpoints = edges(1:end-1) + delta_bin/2;
            
        tickfontsize = 15;
            
        fig = figure;
        hold on;
        
        lgd = gobjects(nTypes,1);
        for iType = 1:nTypes 
            type = cell_types{iType};
            Color = colors{iType};
            switch type
                case 'stellate'
                    shortCellName = 'SC';
                case 'pyramidal'
                    shortCellName = 'PC';
                case 'fast spiking'
                    shortCellName = 'FS';
            end
    
            data = allData{iType};
    
            if ~isempty(data)
                % plot SEM
                meanData = mean(data,1);
                SEM = std(data,1) ./ sqrt(size(data,1));
                meanData = reshape(meanData,1,[]); % column vector
                SEM = reshape(SEM,1,[]); % column vector
                xSEM = [bin_midpoints fliplr(bin_midpoints)] ;         
                ySEM = [meanData+SEM fliplr(meanData-SEM)];
    
                han1 = fill(xSEM,ySEM,Color);
                han1.FaceColor = Color;    
                han1.FaceAlpha = 0.4;      
                han1.EdgeColor = 'none'; 
                drawnow;
    
                name = sprintf('%s',shortCellName);
                % plot mean
                lgd(iType) = plot(bin_midpoints,meanData,'Color',Color,'DisplayName',name);
            end
    
        end

        hold off
        xlim([min(edges) max(edges)])
        ylim([0 0.3]);
        xlabel('Interspike Firing Rate (Hz)')
        ylabel('Counts (1 / Theta Cycle)')
    
        ax = gca;
        ax.YAxis.FontSize = tickfontsize;
        ax.YAxis.FontWeight = 'bold';
        ax.XAxis.FontSize = tickfontsize;
        ax.XAxis.FontWeight = 'bold';

        legend(lgd,'Location','northeastoutside');

    end

end

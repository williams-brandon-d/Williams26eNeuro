function plotThetaPhaseHistDNQX(dataSets, comment_group1,comment_group2, savenum)

nbins = 30;
edges = linspace(-pi,pi,nbins+1);


cell_types = {'stellate','pyramidal','fast spiking'};

% find ID parameters for analysis
params.locations = 'all';
params.experiments = {'currentclamp'};
params.cell_nums = 'all';
params.protocols = {'theta'};

if strcmp(dataSets,'all'); dataSets = {'Camk2','Thy1','PV Transgenic','PV Viral'}; end

nCommentGroups = 2;
nCellTypes = length(cell_types);

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

allData = cell(nCellTypes,nCommentGroups);

for iCell = 1:nCellTypes
    cellType = cell_types{iCell};
    params.cell_types = cell_types(iCell);

    for iComment = 1:nCommentGroups
        if iComment == 1
            params.comments = comment_group1;
        elseif iComment == 2
            params.comments = comment_group2;
        end
        comments = cell2mat(params.comments);

        IDs = getIDs(info,params);
        
        if (isempty(IDs)); fprintf('No %s %s Files Found.',params.cell_types,comments); continue; end
        
        nIDs = numel(IDs);

        HistArray = zeros(nIDs,nbins);

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
            cell_phases = [file.theta_spike_phases{:}]; % nCycles x cycle phase array
            nCycles = numel(file.cycles);
            if ~isempty(cell_phases)
                [HistArray(iID,:),~] = histcounts(cell_phases,edges);
                HistArray(iID,:) = HistArray(iID,:)/nCycles;
            end
            clearvars file;
        end

        figCell = plotHist(HistArray,edges,cellType); % plot figure for each cell type
        pTitle = {sprintf('%s %s %s',dataSet,cellType, comments);'Theta Stim Phase Histogram'};
        sgtitle(pTitle,'Fontweight','bold');
        
        if savenum
            saveFilename = [saveFolder filesep sprintf('%s %s %s Theta Stim Phase Histogram.svg',dataSet,cellType,comments)];
            print(figCell,'-vector','-dsvg',saveFilename);
        end
    
        close all;

        allData{iCell,iComment} = mean(HistArray,1); % save mean cell data
    end

end

figAll = plotAvgHist(allData,edges,cell_types);
pTitle = sprintf('%s Paired Theta Stim Phase Histogram',dataSet);
sgtitle(pTitle,'Fontweight','bold');

if savenum
    saveFilename = [saveFolder filesep sprintf('%s Theta Stim Phase Histogram Average.svg',dataSet)];
    print(figAll,'-vector','-dsvg',saveFilename);
end


end


    function fig = plotHist(HistArray,edges,cellType)

    % plot histogram for each cell type

    delta_bin = edges(2) - edges(1);
    bin_midpoints = edges(1:end-1) + delta_bin/2;
        
    tickfontsize = 15;
    
    switch cellType
        case 'stellate'
            color = [1 0 0];
        case 'pyramidal'
            color = [0 1 0];
        case 'fast spiking'
            color = [0 0 1];
    end
        
    fig = figure;
    hold on;

    plot(bin_midpoints,HistArray,'LineStyle','--','Color',color);
    
    plot(bin_midpoints,mean(HistArray,1),'LineStyle','-','Color',color,'Linewidth',2);
    
    hold off
    ylim([0 1])
    xlim([-pi pi])
    xticks([-pi -pi/2 0 pi/2 pi]);
    xticklabels({'-π','-π/2','0','π/2','π'});
    xlabel('Stim Theta Phase (rad)')
    ylabel('Spike Probability')
    
    ax = gca;
    ax.YTick = 0:0.2:1;
    ax.YAxis.FontSize = tickfontsize;
    ax.YAxis.FontWeight = 'bold';
    ax.XAxis.FontSize = tickfontsize;
    ax.XAxis.FontWeight = 'bold';

end



    function fig = plotAvgHist(allData,edges,cell_types)

    % plot average histograms
    
    nTypes = numel(cell_types);
    nComms = 2;
    
    delta_bin = edges(2) - edges(1);
    bin_midpoints = edges(1:end-1) + delta_bin/2;
    
    tickfontsize = 15;
        
    fig = figure;
    hold on;
    
    for iType = 1:nTypes 
        type = cell_types{iType};

        switch type
            case 'stellate'
                color = [1 0 0];
                shortCellName = 'SC';
            case 'pyramidal'
                color = [0 1 0];
                shortCellName = 'PC';
            case 'fast spiking'
                color = [0 0 1];
                shortCellName = 'FS';
        end
    
        for iComm = 1:nComms
            if iComm == 1
                style = '-';
                shortCommentName = 'CTRL';
            elseif iComm == 2
                style = '--';
                shortCommentName = 'DNQX';
            end

            data = allData{iType,iComm};
    
            if ~isempty(data)
                name = sprintf('%s_{%s}',shortCellName,shortCommentName);
                plot(bin_midpoints,data,'Color',color,'LineStyle',style,'DisplayName',name);
            end
    
        end
    end
    
    hold off
    ylim([0 1])
    xlim([-pi pi])
    xticks([-pi -pi/2 0 pi/2 pi]);
    xticklabels({'-π','-π/2','0','π/2','π'});
    xlabel('Stim Theta Phase (rad)')
    ylabel('Spike Probability')
    
    ax = gca;
    ax.YTick = 0:0.2:1;
    ax.YAxis.FontSize = tickfontsize;
    ax.YAxis.FontWeight = 'bold';
    ax.XAxis.FontSize = tickfontsize;
    ax.XAxis.FontWeight = 'bold';
    
    legend('Location','northeastoutside');

end


end
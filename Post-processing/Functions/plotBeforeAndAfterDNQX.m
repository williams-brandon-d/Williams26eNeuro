function plotBeforeAndAfterDNQX(dataSets, comment_group1,comment_group2, savenum)

fieldName = 'cell';
gammaRange = [50 150]; % Hz for PSD

cell_types = {'stellate','pyramidal','fast spiking'};

% find ID parameters for analysis
params.locations = 'all';
params.experiments = {'inhibition'};
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

        IDs = removeIDs(IDs,info,dataSet); % skip bad recordings - get params for cells to skip
        
        if (isempty(IDs)); fprintf('No %s %s Files Found.',params.cell_types,comments); continue; end
        
        nIDs = numel(IDs);

        zAvg = [];

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
            
            [x,y,z,file.(fieldName).CWTcycleValues] = thetaCWT(file,fieldName, file.cycle_start_index_noArtifacts); % cell average scalogram x,y,z

%             [file.(fieldName).psd,~] = plotCyclePSD(file.(fieldName).raw_data,file.Fs,file.cycle_start_index_noArtifacts,file.cycle_length,file.PSDtype,file.(fieldName).data_units,0);
            [file.(fieldName).psd,~] = plotallCyclePSD(file.(fieldName).raw_data,file.Fs,file.cycle_start_index_noArtifacts,file.cycle_length,file.PSDtype,file.(fieldName).data_units,0);

            psd = mean(file.(fieldName).psd.S,2); % cell average PSD
            freq = file.(fieldName).psd.f;

            if iID == 1
                zAvg = z./nIDs;
                psdAvg = psd./nIDs;
                PSD = zeros(nIDs,length(psd));
            else
                zAvg = z./nIDs + zAvg;
                psdAvg = psd./nIDs + psdAvg;
            end
            PSD(iID,:) = psd;

            clearvars file;
        end

        
        % plot average scalogram
        figCWT = plotScalogram(x,y,zAvg,"",[],'pA','theta');
        clim([0 inf]);
        if iComment == 1
            pTitle = sprintf('%s %s Control Average Scalogram',dataSet,cellType);
            cmax = max(zAvg(:));
        else
            pTitle = sprintf('%s %s DNQX Average Scalogram',dataSet,cellType);
            clim([0 cmax]);
        end
%         pTitle = {sprintf('%s %s %s',dataSet,cellType, comments);'Average Scalogram'};
        sgtitle(pTitle,'Fontweight','bold');
        
        if savenum
            saveFilename = [saveFolder filesep sprintf('%s %s %s Average Scalogram.svg',dataSet,cellType,comments)];
            print(figCWT,'-vector','-dsvg',saveFilename);
        end

        if iComment == 1
            psdControl = psdAvg;
            PSDControl = PSD;
        else
            psdDNQX = psdAvg;
            PSDDNQX = PSD;
        end
    
    end

        % compute and plot PSD from each theta cycle without artifacts
        figPSD = plotPSD(freq,PSDControl,PSDDNQX,gammaRange,cellType);
        pTitle = sprintf('%s %s Average PSD',dataSet,cellType);
        sgtitle(pTitle,'Fontweight','bold');%         % compute and plot PSD from all stim data

        if savenum
            saveFilename = [saveFolder filesep sprintf('%s %s Average PSD DNQX.svg',dataSet,cellType)];
            print(figPSD,'-vector','-dsvg',saveFilename);
        end

end


end


    function fig = plotPSD(freq,Scontrol,Sdnqx,gammaRange,cellType)
    
    switch cellType
        case 'stellate'
            color = [1 0 0];
            semcolor = [1 0.5 0.5];
        case 'pyramidal'
            color = [0 1 0];
            semcolor = [0.5 1 0.5];
        case 'fast spiking'
            color = [0 0 1];
            semcolor = [0.5 0.5 1];
    end
        
    tickfontsize = 15;

    meanScontrol = mean(Scontrol,1);
    meanSdnqx = mean(Sdnqx,1);

    meanScontrol = 10*log10(meanScontrol);
    meanSdnqx = 10*log10(meanSdnqx);

    % plot spectra
    fig = figure;
    hold on;

%     plotSEM(freq,Scontrol,semcolor);
    plot(freq,meanScontrol,'LineWidth',1.5,'Color',color,'LineStyle','-','DisplayName','Control');

%     plotSEM(freq,Sdnqx,semcolor);
    plot(freq,meanSdnqx,'LineWidth',1.5,'Color',[color 0.4],'LineStyle','-','DisplayName','DNQX');

    hold off

    xlabel('Frequency (Hz)');
%     ylabel('PSD (pA^{2}/Hz)');
    ylabel('PSD (dB)');

    xlim(gammaRange);
    box off

    legend('Location','Northeast');

    ax = gca;
    ax.YAxis.FontSize = tickfontsize;
    ax.YAxis.FontWeight = 'bold';
    ax.XAxis.FontSize = tickfontsize;
    ax.XAxis.FontWeight = 'bold';

        function plotSEM(x,data,Color)
            x = reshape(x,1,[]);
            if ~isempty(data)
                % plot Mean +- SEM
                meanData = mean(data,1);
                SEM = std(data,1) ./ sqrt(size(data,1));
    
                xSEM = [x fliplr(x)] ;         
                ySEM = [meanData+SEM fliplr(meanData-SEM)];
    
                han1 = fill(xSEM,ySEM,Color);
                han1.FaceColor = Color;    
                han1.EdgeColor = 'none'; 
                drawnow;
            end
        end

    end


end
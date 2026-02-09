function cyclePPCs(dataSets,dataTypes, savenum)

titlefontsize = 15;
cell_types = {'stellate','pyramidal','fast spiking'};
% cell_types = {'fast spiking'};

% find ID parameters for analysis
params.locations = 'all';
params.experiments = {'currentclamp'};
params.cell_nums = 'all';
params.comments = {'','DNQX before','GBZ before' 'DNQX_GABAzine_before', 'Gabazine before'}; % 'GBZ before' 'DNQX_GABAzine_before', 'Gabazine before', 'DNQX before'
params.protocols = {'theta','theta_2chan'};

if strcmp(dataSets,'all'); dataSets = {'Camk2','Thy1','PV Transgenic','PV Viral'}; end

% load data for figure

nCellTypes = numel(cell_types);
nDataTypes = numel(dataTypes);

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


spike1cell = cell(1,nCellTypes);

for iCell = 1:nCellTypes
    cellType = cell_types{iCell};
    params.cell_types = cell_types(iCell);
    IDs = getIDs(info,params);

    validIDs = removeBadLFP(IDs,info,dataSet); % remove bad LFP recordings (CamK2 only)
    
    if (isempty(IDs)); fprintf('No %s Files Found.',params.cell_types); continue; end
    
    nIDs = numel(IDs);

    d = struct;
    spike1IDcell = cell(nIDs,1); % stim only
    
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
        phases = file.theta_spike_phases; % phases
        FRs = file.FRs;
        cycle_idx_cell = file.theta_spike_indices;

        removeIndices = false(numel(FRs),1);

        removeNum = 2; % num of beginning cycles to remove from analysis
        if removeNum
            removeIndices(1:removeNum) = 1; % add to removeIndices
        end
        
        % remove empty cells
        emptyIdx = cellfun('isempty',phases); % empty idx
        removeIndices = removeIndices | emptyIdx; % add empty cells to removeIndices

        % remove burst spikes for PPC analysis - keep first spike in burst
        FR_threshold = 150; % (Hz) firing rate threshold for burst firing
        for i = 1:numel(phases)
            keepIdx = find(FRs{i} < FR_threshold) + 1;
            keepIdx = reshape(keepIdx,1,[]);
            cycle_phases = phases{i};
            if isempty(cycle_phases); continue; end
            phases{i} = cycle_phases([1 keepIdx]); % keep consecutive spikes < threshold rate
            cycle_idx = cycle_idx_cell{i};
            cycle_idx_cell{i} = cycle_idx([1 keepIdx]);  
        end

        for iType = 1:nDataTypes % stim, theta, gamma
            dataType = dataTypes{iType};

            switch dataType
                case 'stim'
                    d.(dataType).phases{iID} = phases(~removeIndices);

                    % gather average first spike phase from each neuron
%                     idx = ~cellfun('isempty',phases); % non empty idx
%                     out = zeros(size(phases)); % allocate
%                     out(idx) = cellfun(@(v)v(1),phases(idx)); % get first element in non empty cells
%                     out(~out) = []; % remove array values for empty cells

                    out = cellfun(@(v)v(1),phases(~removeIndices)); % get first element in non empty cells
%                     spike1IDcell{iID} = out; % all theta cycles
                    spike1IDcell{iID} = mean(out,'omitnan'); % average non empty per cell

                case {'theta','gamma'}
                     if strcmp(p.protocol,'theta_2chan') % only 2 channel contains LFP
                        % check if bad lfp data
                        if any(strcmp(ID,validIDs))
                            d.(dataType).phases{iID} = getSpikePhase(file,cycle_idx_cell,removeIndices,dataType);
                        else
                            d.(dataType).phases{iID} = {};
                        end
                    else
                        d.(dataType).phases{iID} = {};
                    end
            end

        end
    
        clearvars file;
    end

    for iType = 1:nDataTypes
        dataType = dataTypes{iType};
        dataCell = d.(dataType).phases;

        dataCell(cellfun(@isempty,dataCell)) = []; % remove empty cells

% %         N = size(ppcArray,2);
%         allData = cell(N,1);
%         % setup dataCell using ppcArray
%         for i = 1:N
%             data = ppcArray(:,i);
%             data(isnan(data)) = []; % remove nan
%             allData{i,:} = data;
%         end
% 
%         % spike phase stats - compare spike numbers
%         if nDataTypes > 2
%             comments = cell2mat(params.comments);
            saveFilename = [saveFolder filesep sprintf('%s %s %s Phase Locking stats.xlsx',dataSet,cellType,dataType)];
%             stats = myMultipleIndependentGroupStats(allData,labels,comments,saveFilename,savenum);
% 
%         else 
%             stats.kw = [];
%         end

        % plot phase locking for each cell type
        [figCell,ppcArray,labels] = plotCyclePPCs(dataCell,cellType,saveFilename); % plot figure for each cell type
        titleString = sprintf('%s %s %s Phase Locking',dataSet,cellType,dataType);
        title(titleString,'FontSize',titlefontsize,'Fontweight','bold');
                
        if savenum
            saveFilename = [saveFolder filesep sprintf('%s %s %s Phase Locking.svg',dataSet,cellType,dataType)];
            print(figCell,'-vector','-dsvg',saveFilename);
        end
        
    end

    spike1cell{iCell} = cell2mat(spike1IDcell); % concatenate spike 1 phases for all IDs

end

    % stats for first spike phase between cell types
    if nCellTypes > 2
        comments = cell2mat(params.comments);
        saveFilename = [saveFolder filesep sprintf('%s %s Spike 1 Phase stats.xlsx',dataSet,dataType)];
        myMultipleIndependentGroupStats(spike1cell',cell_types,comments,'linear',saveFilename,savenum);
    end

end


    function phases = getSpikePhase(file,cycle_idx_cell,removeIndices,dataType)
    % load spike indices and lfp phase
%     cycle_idx_cell = file.theta_spike_indices;

    cycle_idx_cell = cycle_idx_cell(~removeIndices);

    switch dataType
        case 'theta'
            lfp_phase = file.lfp.theta_phase;
        case 'gamma' % change bandpass from 60-120 Hz
            gamma_data = filterData(file.lfp.raw_data,file.Fs,[60 120],file.filter_order);
            gamma_phase = angle(hilbert(gamma_data)); % from -pi to pi
            lfp_phase = gamma_phase;
%             lfp_phase = file.lfp.gamma_phase; % 50-200 Hz
    end

    cycle_start_indices = file.cycle_start_index_noArtifacts(~removeIndices(file.cycles_noArtifacts));

    % get spike-lfp phases
    nCycles = numel(cycle_start_indices);
    phases = cell(nCycles,1);

    for icycle = 1:nCycles
        cycle_indices = cycle_start_indices(icycle) + (0:file.cycle_length);
        cycle_phase = lfp_phase(cycle_indices);
        cycle_spike_indices = cycle_idx_cell{icycle,:};
        phases{icycle,:} = cycle_phase(cycle_spike_indices)';
    end

end


end
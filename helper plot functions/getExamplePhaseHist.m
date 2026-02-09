function HistArray = getExamplePhaseHist(ID,info,data_path,edges)

    ID_index = find_index(info,'ID',ID);
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
%     cell_phases = [file.theta_spike_phases{:}]; % nCycles x cycle phase array
%     nCycles = numel(file.cycles);

    phases = file.theta_spike_phases;

    removeNum = 2; % num of beginning cycles to remove from analysis
    removeIndices = false(numel(phases),1);
    if removeNum
        removeIndices(1:removeNum) = 1; % add to removeIndices
    end
    phases = phases(~removeIndices);

    nCycles = numel(phases);

    cell_phases = [phases{:}]; % nCycles x cycle phase array


    if ~isempty(cell_phases)
        [HistArray,~] = histcounts(cell_phases,edges);
        HistArray = HistArray/nCycles;
    end



end
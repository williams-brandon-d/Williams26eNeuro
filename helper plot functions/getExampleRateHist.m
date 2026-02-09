function HistArray = getExampleRateHist(ID,info,data_path,edges)

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
%     cell_FRs = cell2mat(file.FRs); % nCycles x cycle phase array
%     nCycles = numel(file.cycles);

    FRs = file.FRs;
%         cell_FRs = cell2mat(FRs); % nCycles x cycle phase array

    removeNum = 2; % num of beginning cycles to remove from analysis
    removeIndices = false(numel(FRs),1);
    if removeNum
        removeIndices(1:removeNum) = 1; % add to removeIndices
    end
    FRs = FRs(~removeIndices);

    nCycles = numel(FRs);

    cell_FRs = cell2mat(FRs); % nCycles x cycle phase array
%     cell_FRs = [FRs{:}]; % nCycles x cycle phase array

    if ~isempty(cell_FRs)
        [HistArray,~] = histcounts(cell_FRs,edges);
        HistArray = HistArray/nCycles;
    end



end
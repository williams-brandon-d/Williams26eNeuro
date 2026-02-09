clearvars; clc;
cd('C:\Users\brndn\OneDrive\Desktop\White Lab\Matlab Files\');
addpath(genpath('./')); % add all folders and subfolders in cd to path

% process all data first before running this program

% for currentclamp + LFP recordings could choose recording with greatest
% lfp gamma ?

dataSets = {'PV Transgenic'};
savenum = 1;

for iSet = 1:numel(dataSets)
    dataSet = dataSets{iSet};
    
    [info,info_fullname,data_path] = getInfo(dataSet);
    info = removeDuplicates(info,data_path);
    if savenum
        save(info_fullname,'info'); % save info
    end

end



function info = removeDuplicates(info,data_path)

% remove duplicates - choose recording with greatest peak power

% for paired recordings - find before recording with greatest power
% then find the after recording with the same led_input

% make sure that all data is processed beforehand so that all values exist
% if all group data doesnt exist, the first index will be kept

% restrict duplicate removal to specific fieldnames
info_fieldnames = {'comments','experiment','protocol'};
a = struct; % % number of info fieldnames should equal a fieldnames
a.comments = {'','Gabazine before'};
a.experiments = {'inhibition','excitation'};
a.protocols = {'theta','theta_2chan'};

a_fieldnames = fieldnames(a);
for i = 1:numel(info_fieldnames)
    tempFieldarray = a.(a_fieldnames{i});
    tempMask = false(size(info,1),1);
    for ii = 1:numel(tempFieldarray)
        tempMask = tempMask | (strcmp({info.(info_fieldnames{i})},tempFieldarray{ii}))';
    end
    if ~exist('mask','var')
        mask = tempMask;
    else
        mask = mask & tempMask;
    end
end
maskIndices = find(mask);
info_mask = info(mask);

% contstruct table with unique characteristics
fnames = {'location','cell_type','protocol','experiment','cell_num'};
info_table = struct2table(info_mask);
newTable = table;
for i = 1:length(fnames)
    fname = fnames{i};
    newTable.(fname) =  info_table.(fname);
end   

% Find the unique rows, along with indices for identifying the duplicates
[uniqueTableRows,~,indexBackFromUnique] = unique(newTable);
disp(head(uniqueTableRows));

% gather peak gamma power for each group of identical recordings
nGroups = size(uniqueTableRows,1);
for i = 1:nGroups
    group_mask = indexBackFromUnique == i;
    groupMaskIndices = find(group_mask);
    IDs = {info_mask(group_mask).ID};
    nIDs = numel(IDs);
    powers = zeros(nIDs,1);
    led_inputs = powers;

    for iID = 1:nIDs
    
        ID = IDs{iID};
        ID_index = find_index(info_mask,'ID',ID);

        p = info_mask(ID_index);
    
        if isempty(p.comments)
            comment = 'No Comment';
        else
            comment = p.comments;
        end
    
        commentFolder = sprintf('%sresults\\%s\\%s\\%s\\%s\\%s\\%s\\%s',data_path,p.location,p.cell_type,p.cell_num,p.experiment,p.protocol,comment);
    
        dataFolder = getIDFolder(commentFolder,ID);
    
        if ~isempty(dataFolder) 
            dataFilename = [dataFolder filesep 'data.mat'];
            load(dataFilename,'file');
            powers(iID) = file.cell.CWTmaxValues(3); % gather relevant data 
            led_inputs(iID) = file.led_input;
            clearvars file;
        end

    end

    % change info comment of non-max values to duplicate
    % if no values exist then first recording is kept
    [~,max_index] = max(powers); 
    groupMaskMaxIndex = groupMaskIndices(max_index);

    duplicate_mask = group_mask;
    duplicate_mask(groupMaskMaxIndex) = 0; % remove max index from group mask
    finalMask = maskIndices(duplicate_mask); % map indices back to original info struct
    [info(finalMask).comments] = deal(append(p.comments,'Duplicate'));

    % if recording is a before part of a pair, keep after recording with same led input
    if contains(p.comments,'before')
        % find after recording group info indices
        pAfter = p;
        pAfter.comments = replace(p.comments,'before','after');
        [afterIDs,afterID_indices] = find_IDs(info,pAfter);
        % label after recordings with different led_input values as duplicate
        before_led_input = led_inputs(max_index);

        if (before_led_input == 0) || (numel(afterIDs) > 1) 
            % if led_inputs dont exist then choose first after recording to keep
            duplicateMask = true(numel(afterID_indices),1);
            duplicateMask(1) = 0;
        else
            % load after recording led inputs if before led input exists
            after_led_inputs = zeros(numel(afterIDs),1);
            for iAfter = 1:numel(afterIDs)
                afterID = afterIDs{iAfter};
                commentFolderAfter = sprintf('%sresults\\%s\\%s\\%s\\%s\\%s\\%s\\%s',data_path,p.location,p.cell_type,p.cell_num,p.experiment,p.protocol,pAfter.comments);
                dataFolderAfter = getIDFolder(commentFolderAfter,afterID);
                if ~isempty(dataFolderAfter) 
                    dataFilenameAfter = [dataFolderAfter filesep 'data.mat'];
                    load(dataFilenameAfter,'file');
                    after_led_inputs(iAfter) = file.led_input;
                    clearvars file;
                end
            end

            sameLEDmask = after_led_inputs == before_led_input;
            if any(sameLEDmask) % assign unmatched led inputs as duplicates
                duplicateMask = ~sameLEDmask;
            else % if no led inputs match keep first recording
                duplicateMask = true(numel(afterID_indices),1);
                duplicateMask(1) = 0;
            end

        end
        [info(afterID_indices(duplicateMask)).comments] = deal(append(pAfter.comments,'Duplicate'));
    end

end

end
function fullIDfolder = getIDFolder(commentFolder,ID)

% find ID in comment folder - needed if more than one ID with the same params
% could choose file with largest values if more than one ID 
list = dir(commentFolder);
tbl = struct2table(list);
folders = tbl(tbl.isdir,:); % folders in list
if ~isempty(folders)
    namedFolders = folders(~matches(folders.name,[".",".."]),:); % real folders with names
    idFolder = namedFolders(contains(namedFolders.name,[ID ' ']),:).name; % folder name containing ID + a space after
    idFolder = idFolder{1}; % char array not cell
    fullIDfolder = fullfile(commentFolder,idFolder);
else
    fullIDfolder = [];
end

end
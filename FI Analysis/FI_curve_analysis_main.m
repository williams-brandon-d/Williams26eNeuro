%% FI curve analysis

clear variables; close all; clc;
cd('C:\Users\brndn\OneDrive\Desktop\White Lab\Matlab Files\');
addpath(genpath('./')); % add all folders and subfolders in cd to path

savenum = 1; % 0 = don't save, 1 = save info
printnum = 0;

dataSets = {'Thy1'};

% compare mean ISI ratio from all current levels

params.protocols = {'FI'};
params.experiments = {'currentclamp'};

params.locations = 'all';
params.cell_nums = 'all';
 
params.cell_types = {'fast spiking'};
params.comments = {'','DNQX before','DNQX before 10 uM','Gabazine before'}; % combine comment data

tic;
for iSet = 1:numel(dataSets)
dataSet = dataSets{iSet};

[info,info_fullname,data_path] = getInfo(dataSet);

IDs = getIDs(info,params);

if (isempty(IDs)); disp('No Files Found.'); return; end

nIDs = length(IDs);

for iID = 35:nIDs
    ID = IDs{iID};
    ID_index = find_index(info,'ID',ID);
    filename = sprintf('%s.abf',ID);
    fprintf('Analyzing %s,File %d/%d\n',filename,iID,nIDs)
    p = info(ID_index);
    p.dataSet = dataSet;

    % folder save path for results --> location,cell_type,cell_num,experiment,comment,protocol
    if isempty(p.comments)
        comment = 'No Comment';
    else
        comment = p.comments;
    end

    saveFolder = sprintf('%sresults\\%s\\%s\\%s\\%s\\%s\\%s\\%s',data_path,p.location,p.cell_type,p.cell_num,p.experiment,p.protocol,comment,ID);

    fullDataName = fullfile(data_path,filename);
    
    file = analyzeFI(fullDataName,printnum,savenum,saveFolder,p); % analyze FI data 

    if savenum
        file.saveFilename = [file.saveFolder filesep 'data.mat'];
        save(file.saveFilename,'file','-mat','-nocompression');
    end


end

end
toc;

% voltage imaging analysis post processing
% combine data from all FOVs
% script for comparing different theta frequencies across all FOVs

clear variables; close all; clc;
cd('D:/Matlab Files/'); % choose your current directory
addpath(genpath('./')); % add all folders and subfolders in cd to path

savenum = 1;

dataPath = 'D:/TICO Voltage Imaging Data/CaMK2-ChR2-Voltron Data/';

savePath = [dataPath 'Summary Figures/'];

sig = 'symbol';

% theta stimulation frequencies for analysis
% stim_freqs = {'4 Hz','8 Hz','12 Hz','16 Hz'};
stim_freqs = {'8 Hz'};
nFreqs = numel(stim_freqs);

dataFreq = cell(nFreqs,1);

for iFreq = 1:nFreqs
    stim_freq = stim_freqs{iFreq};
    
    % find all data folders that contain stim protocol
    foldernames = getFoldernames(dataPath,stim_freq);
    nFolders = numel(foldernames);
    fprintf('Number of Subfolders Found: %d\n',nFolders);
    
    data = cell(nFolders,nFreqs);
    
    for iFolder = 1:nFolders
        fpath = foldernames{iFolder};
        fprintf('Folder %d: %s\n',iFolder,fpath);
            
        % load saved data
        saveFilename = [fpath filesep 'Processed Data.mat'];
        load(saveFilename,'spikes','lfp','finalTraceIdx');
    
        % if data not found skip analysis
        if ~exist('spikes','var')
            disp('File does not contain spike data.');
            continue;
        end
        if ~exist('lfp','var')
            disp('File does not contain lfp data.');
            continue;
        end
        
        spikes = spikes(finalTraceIdx);
        nNeurons = numel(spikes);
        rate = zeros(nNeurons,1);

        for iNeuron = 1:nNeurons
            mask = spikes(iNeuron).masks;
            nSpikes = sum(mask(lfp.cycle_start_indices_ds(2):lfp.stim_indices_ds(end))); % skip first stim cycle
            nCycles = numel(lfp.cycle_start_indices_ds) - 1;
            rate(iNeuron) = nSpikes/nCycles;
        end

        data{iFolder,iFreq} = rate;

    end
    
    % linearize data for each theta frequency
    dataFreq{iFreq} = cell2mat(data(:,iFreq)); % concatenate column vectors

end

% stats
% if nFreqs > 2
    saveFilename = [savePath 'Average Firing Rate stats.xlsx'];
    stats = myMultipleIndependentGroupStats(dataFreq,stim_freqs,'',saveFilename,savenum);
% else 
%     stats.kw = [];
% end

stats.kw.sig = sig;

% setup y axis for violin plots
y.min = 0;
y.max = 8;
y.dy = 2;
y.ticks = y.min:y.dy:y.max;
y.tickLabels = string(y.ticks);
y.labelstring = 'Firing Rate (Spikes/Cycle)';
y.scale = 'linear';

% get colors for different stim frequencies
% colors = lines(nFreqs);
colors = colormap('lines');
colors = colors(1:nFreqs,:);

alphas = ones(nFreqs,1);

% plot all theta frequencies for comparison - boxplot/violin plots
% fig = plotViolin(rFreq,stim_freqs,y,colors,alphas,stats);
fig = myBoxplot(dataFreq,stim_freqs,y,colors,alphas,struct);

% setup stats in myBoxplot for 4 groups

if savenum
    saveFilename = [savePath filesep 'Average Firing Rate Boxplot.svg'];
    print(fig,'-vector','-dsvg',saveFilename);
end

function foldernames = getFoldernames(rootPath,pattern)
% find all subfolders that have .mat files excluding Processed Data folders
% could instead find all folders that have both .raw and .abf file and Processed Data.mat
filelist = dir(fullfile(rootPath, '**\*Processed Data.mat')); 

% find unique folders
foldernames = unique({filelist.folder})';

% only folders with pattern
foldernames = foldernames(contains(foldernames,pattern));

% if foldername contains pulse - remove
foldernames(contains(foldernames,'pulse','IgnoreCase',true)) = [];

% if foldername contains Bad data - remove
foldernames(contains(foldernames,'Bad data','IgnoreCase',true)) = [];
% if foldername contains lfp only - remove
foldernames(contains(foldernames,'LFP only','IgnoreCase',true)) = [];

end

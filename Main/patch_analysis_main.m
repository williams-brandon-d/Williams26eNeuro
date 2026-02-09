% patch analysis main

% separate script for FI curve analysis for cell identification first?

% fast spiking excitation cell 5

clear variables; close all; clc;
cd('C:\Users\brndn\OneDrive\Desktop\White Lab\Matlab Files\');
addpath(genpath('./')); % add all folders and subfolders in cd to path

savenum = 1; % 0 = don't save, 1 = save info

% change gamma bandpass for PV dataset?
dataSets = {'PV Transgenic'};

% find ID parameters for analysis
params.locations = 'all';
params.cell_types = {'pyramidal'};
params.experiments = {'inhibition'};
params.cell_nums = '8';
% params.comments = {'','Gabazine before','DNQX before','DNQX before 10 uM','GBZ before' 'DNQX_GABAzine_before'}; % 'GBZ before' 'DNQX_GABAzine_before' 'Gabazine before'
params.comments = 'all'; % 'GBZ before' 'DNQX_GABAzine_before'
params.protocols = {'pulse','pulse_2chan'};

if strcmp(dataSets,'all'); dataSets = {'Camk2','Thy1','PV Transgenic','PV Viral'}; end

for iSet = 1:numel(dataSets)
dataSet = dataSets{iSet};

[info,info_fullname,data_path] = getInfo(dataSet);

% info = info(1:204);

IDs = getIDs(info,params);

if (isempty(IDs)); disp('No Files Found.'); return; end

nIDs = numel(IDs);

for iID = 1:nIDs

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
    
    switch p.protocol
        case {'theta','theta_2chan','theta_pulse'} % process sinusoidal stim data
            switch p.experiment
                case {'excitation','inhibition'}
                    % analyze inhibition and excitation currents (plotData,thetaCWT,cyclePSD)
                    % analyze LFP (plotData,thetaCWT,cyclePSD)
                    % current + LFP analysis (cross correlation)
                    file = analyzePatchLFP(fullDataName,1,1,saveFolder,p,dataSet); % analyze LFP

                case 'currentclamp'  % analyze current clamp data
                    % analyze LFP (plotData,thetaCWT,cyclePSD,getPhase,segmentGamma,plotHeatMap,plotGammaBins)
                    % analyze cell  - detectSpikes,rasterPlot,histograms
                    file = analyzeThetaVoltageLFP(fullDataName,1,1,saveFolder,p,dataSet);
                    % spike-stim and spike-LFP analysis, same as voltage imaging
            end
        case {'pulse','pulse_2chan'} % process square stim data
            switch p.experiment
                case {'excitation','inhibition'} % inhibition and excitation currents analysis
                    % analyze inhibition and excitation currents (plotData,pulseCWT,cyclePSD,peak amplitude,sum (integrate) current during stim)
                    % write function for calculating cycleAmplitude (Peak and Sum)?
                    % analyze LFP (plotData,pulseCWT,cyclePSD)    
                    file = analyzePulseCurrentsLFP(fullDataName,1,1,saveFolder,p);
                case 'currentclamp'  % analyze current clamp data
    
                    % spike analysis
            end
    
    
    end

    if savenum
        file.saveFilename = [file.saveFolder filesep 'data.mat'];
        save(file.saveFilename,'file','-mat','-nocompression');
    end

    close all;

end

end

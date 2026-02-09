% patch post processing comparisons

clear variables; close all; clc;
cd('C:\Users\brndn\OneDrive\Desktop\White Lab\Matlab Files\');
addpath(genpath('./')); % add all folders and subfolders in cd to path

dataSets = {'Thy1'};
savenum = 1; % 0 = don't save, 1 = save info

% theta 

% gather stimulation intensities used for voltage and current clamp recordings
% boxplotStimData(dataSets, {'','DNQX before','DNQX before 10 uM','GBZ before','DNQX_GABAzine_before', 'Gabazine before'}, savenum);

% plotLightvsCell() % plot Stimulation intensity (light power) vs gamma current power and frequency

% plotLightvsCellPaired(); % plot paired recordings (power/freq) with different light intensities during PV pulse stim

% Thy1_DorsalvsVentral()

% voltage clamp

% before and after dnqx boxplots
% boxplotBeforeAndAfterDNQX({'Thy1'},{'DNQX before','DNQX before 10 uM'},{'DNQX after','DNQX after 10 uM'}, savenum);

plotBeforeAndAfterDNQX({'Thy1'},{'DNQX before','DNQX before 10 uM'},{'DNQX after','DNQX after 10 uM'}, savenum);

% compare scalogram values
% boxplotCellData(dataSets, {'','DNQX before','DNQX before 10 uM','GBZ before','DNQX_GABAzine_before', 'Gabazine before'}, savenum);

% compare PSD values - harmonics make it difficult to estimate peak frequency in gamma range
% plotCellPSDdata(dataSets, {'','DNQX before','GBZ before','DNQX before 10 uM','DNQX_GABAzine_before', 'Gabazine before'}, savenum);

% getLFPGammaPhase(dataSets, {'','DNQX before','GBZ before','DNQX before 10 uM','DNQX_GABAzine_before', 'Gabazine before'});

% plot theta cwt values over each cycle - synaptic depression
% plotThetaCycleValues(dataSets, {'','DNQX before','DNQX before 10 uM','GBZ before','DNQX_GABAzine_before', 'Gabazine before'}, savenum);

% plot theta amplitude across cycles 
% plotThetaAmplitudeValues(dataSets, {'','DNQX before','DNQX before 10 uM','GBZ before','DNQX_GABAzine_before', 'Gabazine before'}, savenum);

% voltage clamp + lfp

% plotCorrCellData(dataSets, {'','DNQX before','DNQX before 10 uM','GBZ before','DNQX_GABAzine_before', 'Gabazine before'}, savenum); % correlation and lag population comparison

% lfp vs cell data
% plotCamk2inhibitionVsLFPpower(); % plot inhibition vs LFP power

% plotCamk2inhibitionVsLFPfreq(); % plot inhibition vs LFP freq

% plotCamk2inhibitionScalogramVsPSDfreq(); % plot inhibition scalogram freq vs psd

% currentclamp

% plotThetaStimSpikePhaseHist(dataSets, savenum); % histogram spike phases for each cell type

% plotThetaStimSpikeRateHist(dataSets, savenum); % histogram all FRs for each cell type

% plotAverageThetaSpikes(dataSets, savenum); % plot violin with average spikes per theta cycle for each cell type

% plotInterspikeFiringRates(dataSets, savenum); % plot interspike firing rates?

% currentclamp data before and after dnqx has low N values - skip plotting firing rates
% plotThetaPhaseHistDNQX({'Thy1'},{'DNQX before','DNQX before 10 uM'},{'DNQX after','DNQX after 10 uM'}, savenum); % plot paired before and after DNQX comparison for thy1

% current clamp - lfp data

% cyclePPCs(dataSets,{'stim','theta','gamma'}, savenum); % phase locking to stim, theta and gamma theta for camk2 data sets

% cycleVSs(dataSets,{'stim','theta','gamma'}, savenum); % phase locking to stim, theta and gamma theta for camk2 data sets

% plotPatchThetaSpikesPolar(dataSets,{'stim','theta','gamma'}, savenum); % polar plots for spike times relative to theta and gamma


%     % plotSTAs(spikes, lfp, frame_rate, 1); % plot STAs for each neuron for both theta and gamma
%     plotCycleSTAs(spikes, lfp, frame_rate, 1) % compute STA for each first, second, third spike etc.

%     % cycle by cycle spike coherence?
%     plotSpikeCoherence(spikes, lfp, frame_rate, 1); % spike-field coherence

%     spikes = spikeGammaBinHistogram(spikes,lfp,1); % combine cell type data for patch
%     spikeGammaBinHeatMap(spikes,lfp, 1); % spike gamma bin heatmap instead of histogram - sort by cycle #


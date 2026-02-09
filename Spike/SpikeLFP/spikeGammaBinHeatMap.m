function spikeGammaBinHeatMap(spikes,lfp, savenum)
% plot spike probability for binned gamma cycles as heat map for each theta stim cycle 

nSelected = numel(spikes(1).locs);
nGammaBins = numel(spikes.gamma_hist{1});

tickfontsize = 20;
xvalues = spikes.gamma_bin_midpoints;
yvalues = 1:nSelected;

data = zeros(nSelected,nGammaBins);

for iNeuron = 1:nSelected
    data(iNeuron,:) = spikes.gamma_hist{iNeuron};
end

% % skip gamma cycle 0 column
% data = data(:,2:end);
% nGammaBins = size(data,2);
% xvalues = xvalues(2:end);

% setup sorting params
% columns = 1:nGammaBins; % sort by cycle 0 first
columns = 2:nGammaBins; % sort by cycle 1 first
nCols = numel(columns);
direction = 'descend';
directions = cell(1,nCols); % sort all columns by the same direction
[directions{1:nCols}] = deal(direction);

% sort by cluster with highest spike probability
data = sortrows(data,columns,directions);

data = [data; mean(data,1)]; % add mean to bottom of heatmap
yvalues = [yvalues, nSelected+1]; % number of neurons + mean

% plot heatmap of gamma cycle spike probability for all neurons
factor = 10; % plot tick marks every factor
yTicks = factor*(1:1:floor(nSelected/factor)); % tick marks

fig = figure;
imagesc(xvalues,yvalues,data);
% xticks([-pi -pi/2 0 pi/2 pi]);
% xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
yticks([yTicks, size(data,1)]);
yticklabels([string(yTicks), "Avg"]);
colormap('hot');

xlabel('Gamma Cycle #');
ylabel('Neuron #');
title('Spike Rate (1/Gamma Bin/Theta Cycle)','Fontsize',tickfontsize,'Fontweight','bold');
colorbar('Fontsize',tickfontsize,'Fontweight','bold');
clim([0 1]);
% title(h,'','Fontsize',tickfontsize,'Fontweight','bold');

ax = gca;
ax.YAxis.FontSize = tickfontsize;
ax.YAxis.FontWeight = 'bold';
ax.XAxis.FontSize = tickfontsize;
ax.XAxis.FontWeight = 'bold';

if savenum
    saveas(fig,[lfp.dataPath filesep 'Gamma Bin Heatmap.png']);
end

end
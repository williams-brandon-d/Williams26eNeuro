function spikes = spikeGamma(spikes,lfp, savenum)
% for patch combine data from each cell type

nSelected = numel(spikes(1).locs);
nCycles = numel(file.cycles);

% nSamples = numel(lfp.stim_indices_ds);
% phase = linspace(-pi,pi,nSamples);

spikes.gamma_cycle = cell(nSelected,1);
spikes.gamma_hist = cell(nSelected,1);

for iNeuron = 1:nSelected 
    spike_mask = spikes.masks{iNeuron};
    spike_mask_stim = spike_mask(lfp.stim_indices_ds); % restrict time axis to stimulation period 
    spikes.gamma_cycle{iNeuron} = lfp.gamma_cycle_bins(spike_mask_stim);
end

% plot histogram of average
maxbins = max(cellfun(@max,spikes.gamma_cycle));
minbins = min(cellfun(@min,spikes.gamma_cycle));
edges = minbins-0.5:maxbins+0.5;
bin_midpoints = minbins:maxbins;
spikes.gamma_bin_midpoints = bin_midpoints;

% labelfontsize = 20;
tickfontsize = 12;
titlefontsize = 15;

figure;
hold on;

for i = 1:nSelected 

    spikes_gamma = spikes.gamma_cycle{i};

    [N,~] = histcounts(spikes_gamma,edges);
    spike_prob = N/nCycles;
    plot(bin_midpoints,spike_prob,'LineStyle','--','Color',0.5*[1 1 1]);
    
    if i == 1
        all_spikes = spike_prob;
    else

        all_spikes = all_spikes + spike_prob;
    end
    spikes.gamma_hist{i} = spike_prob;
end

plot(bin_midpoints,all_spikes/nSelected,'LineStyle','-','Color','k','Linewidth',2);
hold off

xlim([minbins-0.5 maxbins+0.5])
ylim([-inf inf]);
% ylim([0 1])

% xticks([-pi -pi/2 0 pi/2 pi]);
% xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
xlabel('Gamma Cycle Bins')
ylabel('Spike Rate') % (1/Gamma Bin/Theta Cycle)

% plotTitle = sprintf('%s (n = %d)',lfp.dataPath,nSelected);
plotTitle = sprintf('(n = %d)',nSelected);
sgtitle(plotTitle,'FontSize',titlefontsize,'FontWeight','bold','Interpreter','none');

ax = gca;
% ax.YTick = 0:0.2:1;
ax.YAxis.FontSize = tickfontsize;
ax.YAxis.FontWeight = 'bold';
ax.XAxis.FontSize = tickfontsize;
ax.XAxis.FontWeight = 'bold';
%     ax.YTickLabel = round([0 0.5 1]);

if savenum == 1
    saveName = 'Spike-Gamma Histogram';
    saveas(gcf,[lfp.dataPath filesep saveName '.png']);
end

end
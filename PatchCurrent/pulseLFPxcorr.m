function [file,fig] = pulseLFPxcorr(file)

data_indices = file.pulse_start_index:file.pulse_stop_index;

% get same no artifact trials for both cell and lfp
trials = intersect(file.cell.trials_noArtifacts,file.lfp.trials_noArtifacts);

nSamples = numel(data_indices);
nTrials = numel(trials);

trialCorr = zeros(2*(nSamples-1)+1,nTrials);
lags = trialCorr;

for iTrial = 1:nTrials
    cell_data = file.cell.gamma_data(data_indices,iTrial);
    lfp_data = file.lfp.gamma_data(data_indices,iTrial);
    [trialCorr(:,iTrial), lags] = xcorr(cell_data,lfp_data,'normalized');
end

lags_ms = lags*file.dt*1000; % time in msec
meanCycleCorr = mean(trialCorr,2);
[max_meanCorr,max_meanCorr_index] = max(meanCycleCorr,[],'ComparisonMethod','auto');
max_meanLag = lags_ms(max_meanCorr_index);

[max_peakCorr,max_peakCorr_index] = max(trialCorr(:));
[maxRow,maxCol] = ind2sub(size(trialCorr),max_peakCorr_index);
peakCorr_Lag = lags_ms(maxRow);
maxCorr = trialCorr(:,maxCol);

fontsize = loadFontSizes();
fig = figure;
plot(lags_ms,trialCorr,'Color',0.5*[1 1 1])
hold on
plot(lags_ms,meanCycleCorr,'k','Linewidth',2)
hold off
axis([-50 50 -1 1]);
title(sprintf('Peak Corr. = %.2g    Lag = %.2g ms',max_meanCorr,max_meanLag));
box off
set(gcf, 'Renderer', 'painters');
ax = gca;
ax.XAxis.FontSize = fontsize.tick;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontSize = fontsize.tick;
ax.YAxis.FontWeight = 'bold';   
xlabel('Lag (ms)','FontSize',20,'FontWeight','bold');
ylabel('Correlation Coeff.','FontSize',20,'FontWeight','bold');
title(sprintf('Peak Corr. = %.2g    Lag = %.2g ms',max_meanCorr,max_meanLag),'FontSize',20,'FontWeight','bold');

file.meanLFPcorr = meanCycleCorr;
file.meanLFPcorrCoeff = max_meanCorr;
file.meanLFPcorrLag = max_meanLag;

file.maxLFPcorr = maxCorr;
file.maxLFPcorrCoeff = max_peakCorr;
file.maxLFPcorrLag = peakCorr_Lag;

end
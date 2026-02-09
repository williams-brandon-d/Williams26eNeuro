function [file,fig] = patchLFPxcorr(file)

cycles = 2:6;

nCorrCycles = numel(cycles);

cycleCorr = zeros(2*(file.cycle_length)+1,nCorrCycles);
lags = cycleCorr;

for iCycle = 1:nCorrCycles
    cycle = cycles(iCycle);
    cycle_start = file.cycle_start_index_noArtifacts(cycle);
    cycle_stop = cycle_start + file.cycle_length; 
    cell_data = file.cell.gamma_data(cycle_start:cycle_stop);
    lfp_data = file.lfp.gamma_data(cycle_start:cycle_stop);
    [cycleCorr(:,iCycle), lags] = xcorr(cell_data,lfp_data,'normalized');
end

lags_ms = lags*file.dt*1000; % time in msec
meanCycleCorr = mean(cycleCorr,2);

% restrict xcorr analysis to +- 10 ms 
mask = abs(lags_ms) < 10;
lag_mask = lags_ms(mask);
corr_mask = meanCycleCorr(mask);

switch file.info.experiment
    case 'inhibition'
        [~,max_meanCorr_index] = max(corr_mask,[],'ComparisonMethod','auto'); % find max positive peak
    case 'excitation'
        [~,max_meanCorr_index] = max(-1*corr_mask,[],'ComparisonMethod','auto'); % find max negative peak
end

max_meanLag = lag_mask(max_meanCorr_index);
max_meanCorr = corr_mask(max_meanCorr_index);

[max_peakCorr,max_peakCorr_index] = max(cycleCorr(:));
[maxRow,maxCol] = ind2sub(size(cycleCorr),max_peakCorr_index);
peakCorr_Lag = lags_ms(maxRow);
maxCorr = cycleCorr(:,maxCol);

fontsize = loadFontSizes();
fig = figure;
plot(lags_ms,cycleCorr,'Color',0.5*[1 1 1])
hold on
plot(lags_ms,meanCycleCorr,'k','Linewidth',2)
hold off
axis([-50 50 -1 1]);
box off
set(gcf, 'Renderer', 'painters');
ax = gca;
ax.XAxis.FontSize = fontsize.tick;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontSize = fontsize.tick;
ax.YAxis.FontWeight = 'bold';   
xlabel('Lag (ms)','FontSize',20,'FontWeight','bold');
ylabel('Correlation Coeff.','FontSize',20,'FontWeight','bold');
title(sprintf('Peak Corr. = %.2g    Lag = %.2g ms    ',max_meanCorr,max_meanLag),'FontSize',20,'FontWeight','bold');

file.meanLFPcorr = meanCycleCorr;
file.meanLFPcorrCoeff = max_meanCorr;
file.meanLFPcorrLag = max_meanLag;

file.maxLFPcorr = maxCorr;
file.maxLFPcorrCoeff = max_peakCorr;
file.maxLFPcorrLag = peakCorr_Lag;

end
function fig = spikeGammaOverlay(file)
% plot spike raster on top of gamma bins


delta_phase = 2*pi/file.cycle_length;
cycle_phase = (-pi:delta_phase:pi)'; 

nCycles = numel(file.cycle_start_index_noArtifacts);

% plot params
LineFormat = struct();
switch file.info.cell_type
    case 'stellate'
        LineFormat.Color = [1 0 0];
    case 'fast spiking'
        LineFormat.Color = [0 0 1];  
    case 'pyramidal' 
        LineFormat.Color = [0 1 0];
end
LineFormat.LineWidth = 2;
LineFormat.LineStyle = '-';

% labelfontsize = 20;
tickfontsize = 15;
titlefontsize = 15;




fig = figure;
if any(~emptyCells)
    [~, ~] = plotSpikeRaster(file.cycle_phase_cell,'PlotType','vertline','LineFormat',LineFormat);
end
xlim([cycle_phase(1) cycle_phase(end)])
ylim([0 nCycles]+0.5)
xticks([-pi -pi/2 0 pi/2 pi]);
xticklabels({'-π','-π/2','0','π/2','π'});
xlabel('Stim Theta Phase (rad)')
ylabel('Theta Cycle #')
title(plotTitle,'FontSize',titlefontsize,'FontWeight','bold')
ax = gca;
ax.YAxis.FontSize = tickfontsize;
ax.YAxis.FontWeight = 'bold';
ax.XAxis.FontSize = tickfontsize;
ax.XAxis.FontWeight = 'bold';






end
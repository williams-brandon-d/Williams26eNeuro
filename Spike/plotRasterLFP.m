function fig = plotRasterLFP(file)
% plot average gamma per theta cycle on raster plot

% segment avg cycle lfp gamma 
% draw highlight gamma cycle on plot

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

plotTitle = strjoin({ file.info.cell_type file.info.cell_num sprintf('LED Input: %g mV',file.led_input) });

emptyCells = cellfun('isempty',file.theta_spike_phases);

fig = figure;
sgtitle(plotTitle,'FontSize',titlefontsize,'FontWeight','bold')

ax1 = subplot(2,1,1);
plot(cycle_phase,file.lfp.cycle_gamma_avg,'');
xlim([cycle_phase(1) cycle_phase(end)]);
set(ax1,'visible','off')
drawnow;
% ax1.YAxis.Label.Visible = 'on';
% ylabel(ax1,{'Avg.';'Filtered';'LFP'});
% ax1.YAxis.FontSize = tickfontsize;
% ax1.YAxis.FontWeight = 'bold';
pos1 = get(ax1, 'Position'); % pos array = [x y width height]

ax2 = subplot(2,1,2);
if any(~emptyCells)
    [~, ~] = plotSpikeRaster(file.theta_spike_phases,'PlotType','vertline','LineFormat',LineFormat);
end
xlim([cycle_phase(1) cycle_phase(end)])
ylim([0 nCycles]+0.5)
xticks([-pi -pi/2 0 pi/2 pi]);
xticklabels({'-π','-π/2','0','π/2','π'});
xlabel('Stim Theta Phase (rad)')
ylabel('Theta Cycle #')

ax2.YAxis.FontSize = tickfontsize;
ax2.YAxis.FontWeight = 'bold';
ax2.XAxis.FontSize = tickfontsize;
ax2.XAxis.FontWeight = 'bold';
drawnow;

pos2 = get(ax2, 'Position'); % pos array = [x y width height]

pos1(2) = 0.7;
pos1(4) = 0.2;
set(ax1, 'Position', pos1)

drawnow;

pos2(2) = 0.2;
pos2(4) = 0.5;
set(ax2, 'Position', pos2)


end
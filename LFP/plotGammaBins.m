function plotGammaBins(file, savenum)
% plot lfp gamma cycle as different colors for each stim cycle 

% color = load_fav_colors();

nCycles = numel(file.cycles);

delta_phase = 2*pi/file.cycle_length;
cycle_phase = (-pi:delta_phase:pi)'; 

data = zeros(nCycles,file.cycle_length+1);

for iCycle = 1:nCycles
    gamma_cycle_start = file.lfp.gamma_start_indices{iCycle};
    gamma_cycle_stop = file.lfp.gamma_stop_indices{iCycle};
    for iGamma = 1:numel(gamma_cycle_start)
        gamma_indices = gamma_cycle_start(iGamma):gamma_cycle_stop(iGamma);
        data(iCycle,gamma_indices) = iGamma;
    end
end

tickfontsize = 25;
xvalues = cycle_phase;
yvalues = file.cycles(1):file.cycles(end);

fig = figure;
imagesc(xvalues,yvalues,data);
xticks([-pi -pi/2 0 pi/2 pi]);
xticklabels({'-π','-π/2','0','π/2','π'});

% yticks([1 5 10 15 20]);
% yticklabels({'1','5','10','15','20'});

factor = 5; % tick mark every factor
factorTickMax = floor(file.cycles(end)/factor);
factorTickMin = ceil(file.cycles(1)/factor);
traceTicks = factor*(factorTickMin:factorTickMax);
yticks(traceTicks);
yticklabels(string(traceTicks));

colormap('hot');

xlabel('Stim Theta Phase (rad)');
ylabel('Theta Cycle #');
title('Gamma Cycle Bins','Fontsize',tickfontsize,'Fontweight','bold');
h = colorbar('Fontsize',tickfontsize,'Fontweight','bold');
title(h,{'Gamma';'Cycle #'},'Fontsize',tickfontsize,'Fontweight','bold');

ax = gca;
ax.YAxis.FontSize = tickfontsize;
ax.YAxis.FontWeight = 'bold';
ax.XAxis.FontSize = tickfontsize;
ax.XAxis.FontWeight = 'bold';

if savenum
%     saveas(fig,[file.saveFolder filesep 'LFP Gamma Bins.png']);
    print(fig,'-vector','-dsvg',[file.saveFolder filesep 'LFP Gamma Bins.svg']);
end


end
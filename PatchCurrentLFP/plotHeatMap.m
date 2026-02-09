function file = plotHeatMap(file, fieldName, savenum)
% plot theta and gamma filtered amplitudes as heat map for each theta stim cycle 

switch fieldName
    case 'cell'
        dataTypes = {'Gamma'};
    case 'lfp'
        dataTypes = {'Theta','Gamma'};
end

% dataTypes = {'Theta','Gamma'};
nTypes = numel(dataTypes);

% nCycles = numel(file.cycle_start_index);
nCycles = numel(file.cycles);

delta_phase = 2*pi/file.cycle_length;
cycle_phase = -pi:delta_phase:pi; 
cycle_phase = cycle_phase';

for i = 1:nTypes

    data = zeros(nCycles,file.cycle_length+1);
    type = dataTypes{i};
    
    for icycle = 1:nCycles
        cycle = file.cycles(icycle);
        cycle_indices = file.cycle_start_index(cycle) + (0:file.cycle_length);
        switch type
            case 'Theta'
                data(icycle,:) = file.(fieldName).theta_data(cycle_indices);
            case 'Gamma'
                data(icycle,:) = file.(fieldName).gamma_data(cycle_indices);
        end
    end
    
%     X = [data; mean(data,1)];
    X = data;

    tickfontsize = 20;

    xvalues = cycle_phase;
    yvalues = 1:size(X,1);

    factor = 5;
    yTicks = factor*( 1:1:(floor(nCycles/factor)) );
    
    fig = figure;
    imagesc(xvalues,yvalues,X);
    xticks([-pi -pi/2 0 pi/2 pi]);
    xticklabels({'-π','-π/2','0','π/2','π'});
%     yticks([1 yTicks nCycles+1]);
%     yticklabels(["1" string(yTicks) "Avg"]);
    yticks([1 yTicks]);
    yticklabels(["1" string(yTicks)]);
    colormap('hot');
    
    xlabel('Stim Theta Phase (rad)');
    ylabel('Theta Cycle #');
    title(sprintf('%s %s amplitude',fieldName,type),'Fontsize',tickfontsize,'Fontweight','bold');
    h = colorbar('Fontsize',tickfontsize,'Fontweight','bold');
    title(h,file.(fieldName).data_units,'Fontsize',tickfontsize,'Fontweight','bold');
    
    ax = gca;
    ax.YAxis.FontSize = tickfontsize;
    ax.YAxis.FontWeight = 'bold';
    ax.XAxis.FontSize = tickfontsize;
    ax.XAxis.FontWeight = 'bold';

    if savenum
        print(fig,'-vector','-dsvg',[file.saveFolder filesep sprintf('%s %s Heatmap.svg',fieldName,type)]);
%         saveas(fig,[file.saveFolder filesep sprintf('%s %s Heatmap.svg',fieldName,type)]);
    end

    file.(fieldName).cycle_gamma_avg = mean(data,1);

end

end
function [fig,ppcArray,labels] = plotCyclePPCs(dataCell,cellType,saveFilename)

% types = {'Stim','Theta','Gamma'};

minTotalSpikes = 10;

switch cellType
    case 'stellate'
        maxSpikes = 3;
    case 'pyramidal'
        maxSpikes = 3;
    case 'fast spiking'
        maxSpikes = 4;
end

nNeurons = numel(dataCell);

groupLabel = myGroupLabels({cellType});
groupLabel = groupLabel{1};

[color,~] = getColorAlpha({cellType},1);
color = color{1};

ppc0 = cell(nNeurons,maxSpikes);

for iNeuron = 1:nNeurons
    spike_phases = dataCell{iNeuron}; % nCycles x 1 cell array
    nCycles = numel(dataCell{iNeuron});

     for iSpike = 1:maxSpikes
        phase_cell = cell(nCycles,1); % make raster cell for each trace
        
        for iCycle =  1:nCycles
            cycle_spike_phases = spike_phases{iCycle}; % vector of spike phases

            if numel(cycle_spike_phases) < iSpike
                phase_cell{iCycle} = double.empty(1,0); 
            else
                phase_cell{iCycle} = cycle_spike_phases(iSpike);
            end

        end

        phase_cell = phase_cell(~cellfun('isempty',phase_cell)); % remove empty cells 
        phases = cell2mat(phase_cell);

        % if number of cycles with spikes is < 10... do not record PPC
        if numel(phases) < minTotalSpikes
            ppc0{iNeuron,iSpike} = NaN;
        else
            ppc0{iNeuron,iSpike} = PPC(phases); 
        end

     end
end

ppcArray = cell2mat(ppc0);

% %  boxplots instead? or violin

y.min = 0; y.max = 1; y.dy = 0.2;
y.ticks = y.min:y.dy:y.max;
y.tickLabels = compose('%.1g',y.ticks);
y.labelstring = 'PPC';
y.scale = 'linear';

allData = cell(maxSpikes,1);
Colors = allData;
labels = allData;
alphas = allData;
alpha = 1;
% setup dataCell using ppcArray
for ii = 1:maxSpikes
    data = ppcArray(:,ii);
    data(isnan(data)) = []; % remove nan
    allData{ii,:} = data;
    Colors{ii} = color;
    labels{ii} = sprintf('%s_{Spike%d}',groupLabel,ii);
    alphas{ii} = alpha;
end

% spike phase stats - compare spike numbers
stats = myMultipleIndependentGroupStats(allData,labels,'all','linear',saveFilename,1);
stats.kw.sig = 'symbol'; % symbol or exact

if maxSpikes < 2
    stats.kw = [];
end

% fig = myBoxplot(allData,labels,y,Colors,alphas,[]);
fig = myViolinplot(allData,labels,y,Colors,alphas,stats.kw);
ylim([-0.1 1.1]);

% half violin plots
% [fig,labels] = plotPPC(ppcArray,color,groupLabel); % plot PPC
% ylim([-0.19 1.19]);

    function [fig,labels] = plotPPC(ppc0,color,xlabel)
        % add nice colors
%         c =  [0.45, 0.80, 0.69;...
%               0.98, 0.40, 0.35;...
%               0.55, 0.60, 0.79;...
%               0.90, 0.70, 0.30]; 

        % check if any columns are all NaN - remove group but save correct labels
        nGroups = size(ppc0,2);

        lin_ppc0 = ppc0(:);

        nan_mask = isnan(lin_ppc0);
        lin_ppc0(nan_mask) = []; % remove NaNs from linear array

        groups = ones(numel(ppc0),1);
        for i = 2:nGroups
            idx_start = ((i-1)*size(ppc0,1)) + 1;
            idx_stop = idx_start + size(ppc0,1) - 1;
            groups(idx_start:idx_stop) = i*ones(size(ppc0,1),1);
        end
        groups(nan_mask) = []; % remove NaNs from linear array

        newGroups = unique(groups);
        newNGroups = numel(newGroups);

        colors = repmat(color,[newNGroups, 1]);

        labels = cell(newNGroups,1);
        for iGroup = 1:newNGroups
            group = newGroups(iGroup);
            ppc_group = lin_ppc0(groups == group);
            N = numel(ppc_group);
%             labels{iGroup} = sprintf('%s_{Spike%d} (n=%d)',cellType,group,N);
            if isempty(xlabel)
                labels{iGroup} = sprintf('Spike%d (n=%d)',group,N);
            else
                labels{iGroup} = sprintf('%s_{Spike%d}',xlabel,group);
            end
        end

        tickfontsize = 15;

        % plot PPC
        fig = figure;
        if ~isempty(lin_ppc0)
            if isempty(color)
                daviolinplot(lin_ppc0,'groups',groups,'xtlabels',labels); % default colors
            else
                % vs = violinplot(MPG, Origin);
                daviolinplot(lin_ppc0,'groups',groups,'xtlabels',labels,'color',colors); % my colors
%                 daviolinplot(lin_ppc0,'groups',groups,'xtlabels',labels,'color',colors,'smoothing',0.05); % my colors
            end
        end
        ylim([-0.5 1.5])
        ylabel({'PPC'});
        ax = gca;
        ax.YAxis.FontSize = tickfontsize;
        ax.YAxis.FontWeight = 'bold';
        ax.XAxis.FontSize = tickfontsize;
        ax.XAxis.FontWeight = 'bold';
    end

end
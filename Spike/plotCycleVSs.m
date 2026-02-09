function [figVS,VSArray,figMU,muArray,labels] = plotCycleVSs(dataCell,cellType,dataType,saveFilename)

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

VS = cell(nNeurons,maxSpikes);
mu = VS;

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
            VS{iNeuron,iSpike} = NaN;
            mu{iNeuron,iSpike} = NaN;
        else
            [VS{iNeuron,iSpike},MU,z] = vector_strength(phases); 
            switch dataType
                case {'theta','gamma'}
                    z = -1*z; % shift 180 degrees
                    mu{iNeuron,iSpike} = angle(z);
                    muTickLabels = {'0','π/2','-π/π','-π/2','0'};
                case 'stim'
                    mu{iNeuron,iSpike} = MU;
                    muTickLabels = {'-π','-π/2','0','π/2','π'};
            end
        end

     end
end

VSArray = cell2mat(VS);
muArray = cell2mat(mu);

allVS = cell(maxSpikes,1);
allMU = allVS;

Colors = allVS;
labels = allVS;
alphas = allVS;
alpha = 1;

% setup dataCell using ppcArray
for ii = 1:maxSpikes
    data = VSArray(:,ii);
    data(isnan(data)) = []; % remove nan
    allVS{ii,:} = data;

    mask = data > -inf; % keep angles with VS > threshold
    data = muArray(:,ii);
    data(isnan(data)) = []; % remove nan
    allMU{ii,:} = data(mask);    

    Colors{ii} = color;
    labels{ii} = sprintf('%s_{Spike%d}',groupLabel,ii);
    alphas{ii} = alpha;
end

% spike phase stats - compare spike numbers
VSstats = myMultipleIndependentGroupStats(allVS,labels,'all','linear',saveFilename,1);
VSstats.kw.sig = 'symbol'; % symbol or exact

MUstats = myMultipleIndependentGroupStats(allMU,labels,'all','linear',saveFilename,1);
MUstats.kw.sig = 'symbol'; % symbol or exact

if maxSpikes < 2
    VSstats.kw = [];
    MUstats.kw = [];
end

% plot VS
y.min = 0; y.max = 1; y.dy = 0.2;
y.ticks = y.min:y.dy:y.max;
y.tickLabels = compose('%.1g',y.ticks);
y.labelstring = 'VS';
y.scale = 'linear';

% fig = myBoxplot(allData,labels,y,Colors,alphas,[]);
figVS = myViolinplot(allVS,labels,y,Colors,alphas,VSstats.kw);
ylim([-0.1 1.1]);

% half violin plots
% [fig,labels] = plotPPC(ppcArray,color,groupLabel); % plot PPC
% ylim([-0.19 1.19]);

% plot MU
y.min = -pi; y.max = pi; y.dy = pi/2;
% y.min = 0; y.max = 360; y.dy = 90;
y.ticks = y.min:y.dy:y.max;
y.tickLabels = muTickLabels;
% y.tickLabels = compose('%.1g',y.ticks);
y.labelstring = 'Phase (rad)';
y.scale = 'linear';

% fig = myBoxplot(allData,labels,y,Colors,alphas,[]);
figMU = myViolinplot(allMU,labels,y,Colors,alphas,MUstats.kw);

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
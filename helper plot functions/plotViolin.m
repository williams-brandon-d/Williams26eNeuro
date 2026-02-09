function fig = plotViolin(allData,xTickLabels,y,colors,alphas,stats)
    % add nice colors
%         c =  [0.45, 0.80, 0.69;...
%               0.98, 0.40, 0.35;...
%               0.55, 0.60, 0.79;...
%               0.90, 0.70, 0.30]; 

    % check if any columns are all NaN - remove group but save correct labels
    nGroups = numel(allData);

    % pad with NaN
    maxNumEl = max(cellfun(@numel,allData(:)));
    data_all_pad = cellfun(@(x){padarray(x(:),[maxNumEl-numel(x),0],NaN,'post')}, allData(:)); % Pad each vector with NaN values to equate lengths
    dataArray = cell2mat(data_all_pad'); 

    linData = dataArray(:);

    nan_mask = isnan(linData);
    linData(nan_mask) = []; % remove NaNs from linear array

    groups = ones(numel(dataArray),1);
    for i = 2:nGroups
        idx_start = ((i-1)*size(dataArray,1)) + 1;
        idx_stop = idx_start + size(dataArray,1) - 1;
        groups(idx_start:idx_stop) = i*ones(size(dataArray,1),1);
    end
    groups(nan_mask) = []; % remove NaNs from linear array

    newGroups = unique(groups);

    colors = cell2mat(colors(newGroups));

%     newNGroups = numel(newGroups);
% %         labels = cell(newNGroups,1);
%     row1 = cell(1,newNGroups);
%     row2 = row1;
%     for iGroup = 1:newNGroups
%         group = newGroups(iGroup);
%         data_group = linData(groups == group);
%         N = numel(data_group);
%         row1{iGroup} = sprintf('%s',cell_types{iGroup});
%         row2{iGroup} = sprintf('(n=%d)',N);
%     end
%     labelArray = [row1;row2];
%     xTickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));

    tickfontsize = 15;

    % if yaxis is log scale - take log of data and plot on linear scale
    if strcmp(y.scale,'log')
        linData = log10(linData);
    end

    % plot violin
    fig = figure;
    if ~isempty(linData)
        if nargin < 3 
            daviolinplot(linData,'groups',groups); % default colors
        else
            daviolinplot(linData,'groups',groups,'color',colors,'violinalpha',alphas{1},'scatter',1); % my colors and alphas
        end
    end

    ax = gca;
    ax.TickLabelInterpreter = 'tex';
    ax.XTickLabel = xTickLabels;
    ax.XAxis.FontSize = tickfontsize;
    ax.XAxis.FontWeight = 'bold';  
    ylabel(y.labelstring,'FontSize',25,'FontWeight','bold')
    
    if strcmp(y.scale,'log')
        ax.YLim = log10([y.min y.max]);
        ax.YTick = log10(y.ticks);
    else
        ax.YLim = [y.min y.max];
        ax.YTick = y.ticks;
    end

%     ax.YAxis.Scale = y.scale; 
    ax.YTickLabel = y.tickLabels;
    ax.YAxis.FontSize = tickfontsize;
    ax.YAxis.FontWeight = 'bold';

    xPos = 1:nGroups;

    if isfield(stats,'results')
        if isfield(stats,'tbl')
            pvalue = stats.tbl{1,end};
            pvalue = pvalue{1};
            if pvalue > 0.05 % check pvalue for anova or kw test
                return;
            end
        end
        if strcmp(y.scale,'log')
            dataMax = log10(y.max);
            ax.YLim = [log10(y.min) log10(y.max)+1];
        else
            dataMax = y.max;
            ax.YLim = [y.min 1.2*y.max];
        end
        plotSignificance(stats,xPos,dataMax,ax.YAxis.Scale)
    end


end
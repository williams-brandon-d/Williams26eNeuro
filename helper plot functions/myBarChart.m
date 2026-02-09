function fig = myBarChart(data_all,groupLabels,y,colors,alphas,stats)

tickfontsize = 15;

[nCells,nComms] = size(data_all);

data_trans = reshape(data_all,nComms,[]); % for easier grouping transpose so groups are in order
maxNumEl = max(cellfun(@numel,data_trans(:)));
data_all_pad = cellfun(@(x){padarray(x(:),[maxNumEl-numel(x),0],NaN,'post')}, data_trans(:)); % Pad each vector with NaN values to equate lengths
data_all_mat = cell2mat(data_all_pad'); 

meanData = mean(data_all_mat,1,"omitnan");

N = sum(~isnan(data_all_mat),1);

SEM = std(data_all_mat,1,"omitnan")./sqrt(N);

xPos = 1:nCells;

barWidth = 0.8;

fig = figure;
hold on
if ~isempty(meanData)
    b = bar(xPos,meanData,barWidth,'FaceColor','flat','EdgeColor','none');
    for i = 1:nCells
        b.CData(i,:) = colors{i}; % CData = nCells x RGB matrix
        b.FaceAlpha = alphas{1}; % just one value for group
    end
    er = errorbar(xPos,meanData,SEM);
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
    er.CapSize = 15;
    er.LineWidth = 1;
end
hold off
box off

ax = gca(); 
ax.TickLabelInterpreter = 'tex';  % needed for some plots like boxplot.

xticks(xPos)
xticklabels(groupLabels)
ax.XLim = [xPos(1)-barWidth xPos(end)+barWidth];
ax.XAxis.FontSize = tickfontsize;
ax.XAxis.FontWeight = 'bold';
% xtickangle(ax,45)

ax.YLim = [y.min y.max];

ylabel(y.labelstring,'FontSize',25,'FontWeight','bold')

if contains(y.labelstring,'\Delta')
    switch y.labelstring
        case '\Delta Log Peak Gamma Power (pA^{2})'
            ax.YLim = [-3 1];
        case '\Delta Frequency (Hz)'
            ax.YLim = [-50 50];
        case '\Delta Theta Phase (rad)'
            ax.YLim = [-pi pi];
    end
else
    switch y.labelstring
        case 'Log Peak Gamma Power (pA^{2})'
%             ax.YLim = [1 4]; % thy1 and pv barcharts
            ax.YLim = [1 4]; % camk2 barcharts
        case 'Frequency (Hz)'
            ax.YLim = [50 150];
        case 'Theta Phase (rad)'
            ax.YLim = [-pi pi];
    end

end
ax.YAxis.Scale = y.scale; 
ax.YTick = y.ticks;
ax.YTickLabel = y.tickLabels;
ax.YAxis.FontSize = tickfontsize;
ax.YAxis.FontWeight = 'bold';

if isfield(stats,'results')
    if isfield(stats,'tbl')
        pvalue = stats.tbl{1,end};
        pvalue = pvalue{1};
        if pvalue > 0.05 % check pvalue for anova or kw test
            return;
        end
    end
%     dataMax = y.max;
%     ax.YLim = [y.min 1.2*y.max];
    dataMax = max(meanData+SEM);
    plotSignificance(stats,xPos,dataMax,y.scale)
end

end

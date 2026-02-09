function fig = myAlternativeBoxplot(data_all,xTickLabels,y,colors,alphas,stats)

colors = cell2mat(colors);

tickfontsize = 15;
scattersize = 70;

[nCells,nComms] = size(data_all);

xPos = 1:nCells;

data_trans = data_all'; % for easier grouping transpose so groups are in order
maxNumEl = max(cellfun(@numel,data_trans(:)));
data_all_pad = cellfun(@(x){padarray(x(:),[maxNumEl-numel(x),0],NaN,'post')}, data_trans(:)); % Pad each vector with NaN values to equate lengths
data_all_mat = cell2mat(data_all_pad); 

data_groups = zeros(numel(data_all_mat),1);
count = 0;
for ii = 1:nCells
    for iii = 1:nComms
        count = count + 1;
        ix1 = (count-1).*(maxNumEl) + 1;
        ix2 = ix1 + maxNumEl - 1;
        data_groups(ix1:ix2) = count.*ones(maxNumEl,1);
    end
end


fig = figure;

if ~isempty(data_all_mat)
%     plot_handle = boxplot(data_all_mat,data_groups,'Colors','k','Symbol','k.','OutlierSize',10,'Positions',xPos);

%     h = daboxplot(data_all_mat,'groups',data_groups,'xtlabels', xTickLabels,...
%         'colors',colors,'boxalpha',0.7,'whiskers',0,'scatter',1,'outsymbol','k*',...
%         'outliers',1,'scattersize',25,'flipcolors',0,'boxspacing',1.2,'scattercolors',{'k','w'}); 

    h = daboxplot(data_all_mat,'groups',data_groups,'xtlabels', xTickLabels,...
        'colors',colors,'fill',1,'boxalpha',0.4,'whiskers',0,'scatter',1,'outsymbol','k*',...
        'outliers',1,'scattersize',scattersize,'flipcolors',0,'boxspacing',1.2,'scattercolors',{'k','w'}); 

    for g = 1:numel(xTickLabels)
       set(h.sc(:,g),'MarkerFaceColor',colors(g,:));
       set(h.ot(:,g),'CData',colors(g,:));
    end


end

% box off

% xticks(xPos)
% xticklabels(xTickLabels)

ax = gca(); 
ax.TickLabelInterpreter = 'tex';  % needed for some plots like boxplot.

ax.XAxis.FontSize = tickfontsize;
ax.XAxis.FontWeight = 'bold';
% xtickangle(ax,45)

ylabel(y.labelstring,'FontSize',25,'FontWeight','bold')
ax.YLim = [y.min y.max];
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
        
    dataMax = max(data_all_mat(:));
    plotSignificance(stats,xPos,dataMax,y.scale)
end

end
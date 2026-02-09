function fig = myBoxplot(data_all,xTickLabels,y,colors,alphas,stats)

colors = flip(colors);
alphas = flip(alphas);

tickfontsize = 15;

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
    plot_handle = boxplot(data_all_mat,data_groups,'Colors','k','Symbol','k.','OutlierSize',10,'Positions',xPos);
    set(findobj(gcf,'-regexp','Tag','\w*Whisker'),'LineStyle','-')
    set(plot_handle, 'Linewidth',1.5)
    xticks(xPos)
    xticklabels(xTickLabels)
    h = findobj(gca,'Tag','Box');
    h2 = findobj(gcf, 'type', 'line', 'Tag', 'Median');
    for j=1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),colors{j},'FaceAlpha',alphas{j});
        line(h2(j).XData,h2(j).YData,'Color','k','Linewidth',1); % median 
    end
    drawnow;
%     set(findobj(gcf, 'type', 'line', 'Tag', 'Median'), 'Color', 'k');    
end

box off

% plot data markers
% hold on
% count = 0;
% for ii = 1:nComms
%     for iii = 1:nCells
%           count = count + 1;
%           plot(xPos(count)*ones(length(data_all_pad{count}),1),data_all_pad{count},'xk','Markersize',5)
%     end
% end
% hold off

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
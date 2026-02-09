function fig = myViolinplot(data_all,xTickLabels,y,colors,alphas,stats)

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
%     violinplot(data_all_mat,data_groups,'ViolinColor',cell2mat(colors),'ViolinAlpha',0.5,'BoxColor',0*[1 1 1],'EdgeColor',0*[1 1 1]);
    violinplot(data_all_mat,data_groups,'Bandwidth',0.1,'ViolinColor',cell2mat(colors),'ViolinAlpha',0.5,'BoxColor',0*[1 1 1],'EdgeColor',0*[1 1 1]);
    xticks(xPos)
    xticklabels(xTickLabels)
    drawnow;
end

box off

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
    dataMax = max(data_all_mat(:));
    plotSignificance(stats,xPos,dataMax,y.scale)
end

end
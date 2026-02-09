function t2 = myUnPairedTtestStats(data_all,rowLabels)

[nRows,nCols] = size(data_all);

if (nRows == 1) || (nCols == 1)
    data_all = reshape(data_all,1,[]);
    if nRows == 1
        rowLabels = cell2mat(reshape(rowLabels,1,[]));
    end
else
    rowLabels = reshape(rowLabels,[],1);
end

% [~,t2.p,~,stats] = ttest2(data_all{1,1}, data_all{1,2},'tail','both','Vartype','equal');
% t2.tstat = stats.tstat;
% t2.df = stats.df;

% t2.results = cell2table([cell_types(1) cell_types(2) {t2.tstat} {t2.df} {t2.p}],"VariableNames",{'Group A','Group B','T','df','P-value'});

[nRows,~] = size(data_all);

t2 = struct();
t2.H = cell(nRows,1);
t2.p = t2.H;
t2.tstat = t2.H;
t2.df = t2.H;

% or use cellfun
for i = 1:nRows
    if (numel(data_all{i,1}) > 2) && (numel(data_all{i,2}) > 2)
        [~,t2.p{i},~,stats] = ttest2(data_all{i,1}, data_all{i,2},'tail','both','Vartype','equal');
        t2.tstat{i} = stats.tstat;
        t2.df{i} = stats.df;
    else
        t2.p{i} = NaN; t2.tstat{i} = NaN; t2.df{i} = NaN;
    end
end

t2.results = cell2table([rowLabels t2.p t2.tstat t2.df],'VariableNames',{'Group A vs B','P-value','tstat','df'});

end
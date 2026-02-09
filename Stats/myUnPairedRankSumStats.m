function rs = myUnPairedRankSumStats(data_all,rowLabels)

[nRows,nCols] = size(data_all);

if (nRows == 1) || (nCols == 1)
    data_all = reshape(data_all,1,[]);
    if nRows == 1
        rowLabels = cell2mat(reshape(rowLabels,1,[]));
    end
else
    rowLabels = reshape(rowLabels,[],1);
end

% [rs.p,~,stats] = ranksum(data_all{1,1},data_all{1,2},'method','exact','tail','both');
% rs.w = stats.ranksum;
% 
% rs.results = cell2table([cell_types(1) cell_types(2) {rs.w} {rs.p}],"VariableNames",{'Group A','Group B','W','P-value'});

[nRows,~] = size(data_all);

rs = struct();
rs.p = cell(nRows,1);
rs.w = rs.p;

% or use cellfun
for i = 1:nRows
    if (numel(data_all{i,1}) > 2) && (numel(data_all{i,2}) > 2)
        [rs.p{i},~,stats] = ranksum(data_all{i,1},data_all{i,2},'tail','both');
        rs.w{i} = stats.ranksum;
    else
        rs.p{i} = NaN; rs.w{i} = NaN;
    end
end

% rs.results = cell2table([rowLabels(1) rowLabels(2) rs.p rs.w],'VariableNames',{'Group A','Group B','P-value','w'});
rs.results = cell2table([rowLabels rs.p rs.w],'VariableNames',{'Group A vs B','P-value','w'});

end
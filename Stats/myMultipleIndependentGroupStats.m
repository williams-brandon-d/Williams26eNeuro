function stats = myMultipleIndependentGroupStats(allData,rowLabels,colLabels,scale,saveFilename,savenum)
% compare between independent groups (rowLabels)

nRows = numel(rowLabels);

% log transform data?
if strcmp(scale,'log')
    stats.summary = myLogSummaryTable(allData,rowLabels); % calculate  geometric mean, SEM, and CI intervals
    allData = cellfun(@log10,allData,'UniformOutput',false); % returns cell array
else
    stats.summary = mySummaryTable(allData,rowLabels); % calculate mean, SEM and N 
end

stats.data_table = myDataTable(allData,rowLabels);
stats.sw = mySWstats(allData,rowLabels,colLabels);
stats.levene = myVarTest(allData,rowLabels,colLabels);

[nRows,nCols] = size(allData);
if (nRows == 1) || (nCols == 1)
    nRows = max(nRows,nCols);
end

if nRows > 2
    stats.anova = myAnovaStats(allData,rowLabels);
    stats.kw = myKruskalWallisStats(allData,rowLabels);
elseif nRows == 2 % only 2 groups 
    stats.t2 = myUnPairedTtestStats(allData,rowLabels);
    stats.rs = myUnPairedRankSumStats(allData,rowLabels);
end

if savenum
    writetable(stats.data_table,saveFilename,'Sheet','Data','WriteMode','overwritesheet');  % save data table
    writetable(stats.summary,saveFilename,'Sheet','Summary','WriteMode','overwritesheet','WriteRowNames',true);  % save summary stats table
    writetable(stats.sw.results,saveFilename,'Sheet','Shapiro-Wilk','WriteMode','overwritesheet');  % save sw stats table
    writetable(stats.levene.results,saveFilename,'Sheet','Levene','WriteMode','overwritesheet');  % save levene stats table
    if nRows > 2
        writetable(stats.anova.results,saveFilename,'Sheet','ANOVA-Tukey','WriteMode','overwritesheet');  % save anova stats table
        writetable(stats.anova.tbl,saveFilename,'Sheet','ANOVA','WriteMode','overwritesheet');  % save anova stats table
        writetable(stats.kw.results,saveFilename,'Sheet','Kruskal-Wallis-Dunn-Sidak','WriteMode','overwritesheet');  % save kruskal-wallis stats table
        writetable(stats.kw.tbl,saveFilename,'Sheet','Kruskal-Wallis','WriteMode','overwritesheet');  % save kruskal-wallis stats table
    elseif nRows == 2
        writetable(stats.t2.results,saveFilename,'Sheet','Independent T-test','WriteMode','overwritesheet');  % save independent t-test stats table
        writetable(stats.rs.results,saveFilename,'Sheet','Rank Sum test','WriteMode','overwritesheet');  % save rank sum test stats table
    end
end


end
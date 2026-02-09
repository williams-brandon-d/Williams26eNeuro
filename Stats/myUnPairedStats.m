function stats = myUnPairedStats(allData,cell_types,dataSets,saveFilename,savenum)

stats.data_table1 = myDataTable(allData(:,1),cell_types);
stats.data_table2 = myDataTable(allData(:,2),cell_types);

stats.sw = mySWstats(allData,cell_types,dataSets);

stats.var = myVarTest(allData,cell_types,dataSets);

stats.t2 = myUnPairedTtestStats(allData,cell_types);

stats.rs = myUnPairedRankSumStats(allData,cell_types);

if savenum
    if ~isempty(stats)
        writetable(stats.data_table1,saveFilename,'Sheet',sprintf('%s Data',dataSets{1}),'WriteRowNames',true,'WriteMode','overwritesheet');  % save data table
        writetable(stats.data_table2,saveFilename,'Sheet',sprintf('%s Data',dataSets{2}),'WriteRowNames',true,'WriteMode','overwritesheet');  % save data table
        writetable(stats.t2.results,saveFilename,'Sheet','Independent T-test','WriteRowNames',true,'WriteMode','overwritesheet');  % save t2 stats table
        writetable(stats.rs.results,saveFilename,'Sheet','Rank Sum','WriteRowNames',true,'WriteMode','overwritesheet');  % save sr stats table
        writetable(stats.sw.results,saveFilename,'Sheet','Shapiro-Wilk','WriteRowNames',true,'WriteMode','overwritesheet');  % save sw stats table
        writetable(stats.var.results,saveFilename,'Sheet','Levene','WriteRowNames',true,'WriteMode','overwritesheet');  % save var stats table
    end
end

end
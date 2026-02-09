function a = myAnovaStats(data_all,cell_types)

% construct data array with NaNs
maxNumEl = max(cellfun(@numel,data_all(:)));
data_all_pad = cellfun(@(x){padarray(x(:),[maxNumEl-numel(x),0],NaN,'post')}, data_all(:)); % Pad each vector with NaN values to equate lengths
data_all_mat = cell2mat(data_all_pad'); 

% run one-way anova and multiple comparisons test 
a.testType = 'One-way ANOVA';
[~,tbl,stats] = anova1(data_all_mat,cell_types,'off');

a.tbl = cell2table(tbl(2:end,:),'VariableNames',tbl(1,:));

% if variance is equal (levene's test) - use tukey-kramer
% assuming variance is equal if using anova
% tukey-kramer adjusts for unequal sample sizes
% tukey-kramer is recommended when making many pairwise comparisons
[results,~,~,gnames] = multcompare(stats,CriticalValueType="tukey-kramer",Alpha=0.05,Display="off");

% if variance not equal - use bonferroni or maybe Games-Howell?
% bonferroni is more powerful for a small number of comparisons but usually too
% conservative for large number of comparisons
% bonferroni is recommended when choosing a specific number of comparisons
% rather than making all pairwise comparisons
% [results,~,~,gnames] = multcompare(a.stats,CriticalValueType="bonferroni",Alpha=0.05,Display="off");

% construct table of multiple comparisons results
a.results = array2table(results,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
a.results.("Group A") = gnames(a.results.("Group A"));
a.results.("Group B") = gnames(a.results.("Group B"));

end
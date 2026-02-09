function kw = myKruskalWallisStats(data_all,cell_types)

% construct data array with NaNs
maxNumEl = max(cellfun(@numel,data_all(:)));
data_all_pad = cellfun(@(x){padarray(x(:),[maxNumEl-numel(x),0],NaN,'post')}, data_all(:)); % Pad each vector with NaN values to equate lengths
data_all_mat = cell2mat(data_all_pad'); 

% run one-way anova and multiple comparisons test 
kw.testType = 'One Way Kruskal-Wallis ANOVA';
[~,tbl,stats] = kruskalwallis(data_all_mat,cell_types,'off');

kw.tbl = cell2table(tbl(2:end,:),'VariableNames',tbl(1,:));

% need dunn's test for non-parametric - compute z statistic and apply multiple comparisons adjustment (Dunn-Sidak or Bonferroni)
% bonferroni is most conservative - usually too conservative
% The Bonferroni adjustment multiplies each p-value by m, p∗ = pm
% dunn-sidak adjustment corrects the Bonferroni adjustment's error by defining the family wise error rate and gives a slightly smaller p∗, p∗ = 1−(1−p)m
% p value is only slightly lower for dunn-sidak at low number of comparisons
% [results,~,~,gnames] = multcompare(stats,CriticalValueType="dunn-sidak",Alpha=0.05,Display="off");
[results,~,~,gnames] = multcompare(stats,CriticalValueType="bonferroni",Alpha=0.05,Display="off");

% the math for z statistic and pval is the same for multcompare
% Qtable = myDunnTest(data_all,cell_types); % dunn test with sidak adjustment

% construct table of multiple comparisons results
kw.results = array2table(results,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
kw.results.("Group A") = gnames(kw.results.("Group A"));
kw.results.("Group B") = gnames(kw.results.("Group B"));

%     function Qtable = myDunnTest(data_all,cell_types)
%         x = vertcat(data_all{:});
%         x = reshape(x,1,[]);
%         g = zeros(1,numel(x));
%         Ns = cellfun(@numel,data_all(:));
%         count = 0;
%         for i = 1:numel(cell_types)
%             nGroup = Ns(i);
%             if i == 1
%                 g(1:nGroup) = repmat(i,nGroup,1);
%             else
%                 g(count+(1:nGroup)) = repmat(i,nGroup,1);
%             end
%             count = count + nGroup;
%         end
%         Qtable = dunn(x,g);
% 
%     end

end
function test = myVarTest(data_all,rowLabels,colLabels)

[nRows,nCols] = size(data_all);

if (nRows == 1) || (nCols == 1)
    data_all = reshape(data_all,[],1);
    rowArray = cell2mat(reshape(rowLabels,1,[]));
%     if nRows == 1
%         rowArray = cell2mat(reshape(rowLabels,1,[]));
%     end
else
    rowArray = reshape(rowLabels,[],1);
end

[nRows,nCols] = size(data_all);

% maxSize = max(size(data_all));

test = struct;

test.p = cell(nCols,1);
test.fstat = test.p;
test.df = test.p;

colArray = cell(nCols,1);

% could just pad all cells beforehand
for iCol = 1:nCols
    data_col = data_all(:,iCol);

    % construct data array with NaNs
    maxNumEl = max(cellfun(@numel,data_col(:)));
    data_all_pad = cellfun(@(x){padarray(x(:),[maxNumEl-numel(x),0],NaN,'post')}, data_col(:)); % Pad each vector with NaN values to equate lengths
    data_all_mat = cell2mat(data_all_pad'); % N x Groups
    
    [test.p{iCol},stats] = vartestn(data_all_mat,'TestType','LeveneAbsolute','Display','off');
    test.fstat{iCol} = stats.fstat; 
    test.df{iCol} = stats.df; 

    if nCols == 1
        if iscell(colLabels)
            colArray{iCol} = cell2mat(reshape(colLabels,1,[]));
        else
            colArray{iCol} = colLabels;
        end
    else
        if iscell(colLabels)
            colArray{iCol} = colLabels{iCol};
        else
            colArray{iCol} = colLabels;
        end
    end

end

test.results = cell2table([rowArray colArray test.p test.fstat test.df],"VariableNames",{'Cell Types','Comments','P-value','F','df'});

end
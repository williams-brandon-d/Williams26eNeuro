function data_table = mySummaryTable(data_all,cell_types)
[nCells,nGroups] = size(data_all);

% construct data array with NaNs
maxNumEl = max(cellfun(@numel,data_all(:)));
data_all_pad = cellfun(@(x){padarray(x(:),[maxNumEl-numel(x),0],NaN,'post')}, data_all); % Pad each vector with NaN values to equate lengths
data_all_pad = data_all_pad';

if nGroups > 1
    data_all_pad = data_all_pad(:);
    data = reshape(data_all_pad,1,[]);
    data_all_mat = cell2mat(data); 
else
    data_all_mat = cell2mat(data_all_pad); 
end

% if size(data_all,2) == 1
%     data_all_mat = data_all_mat';
% end

meanData = mean(data_all_mat,1,"omitnan");

N = sum(~isnan(data_all_mat),1);

if all(N < 2)
    sd = NaN*zeros(1,nCells);
else
    sd = std(data_all_mat,1,"omitnan");
end

SEM = sd./sqrt(N);

medianData = median(data_all_mat,1,'omitnan');

Q = prctile(data_all_mat,[25 75],1);

dataMin = min(data_all_mat,[],1,'omitnan');
dataMax = max(data_all_mat,[],1,'omitnan');
dataRange = range(data_all_mat,1);

% column names for cell types
varNames = cell(nCells*nGroups,1);
count = 0;
for i = 1:nCells
    switch cell_types{i}
        case {'stellate','Stellate'}
            cellName = 'Stellate';
        case {'pyramidal','Pyramidal'}
            cellName = 'Pyramidal';
        case {'fast spiking','FastSpiking'}
            cellName = 'FastSpiking';
        otherwise
            cellName = cell_types{i};
    end
    for iGroup = 1:nGroups
        count = count+1;
        groupName = sprintf('group%d',iGroup);
        varNames{count} = [cellName '_' groupName];
    end
end

data_table = array2table([meanData; sd; SEM; N; medianData; Q(1,:); Q(2,:); dataMin; dataMax; dataRange],'VariableNames',varNames,'RowNames',{'Mean','SD','SEM','N','Median','Q25','Q75','Min','Max','Range'});

end
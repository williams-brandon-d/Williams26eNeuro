function data_table = myDataTable(data_all,cell_types)
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

% build data table
data_table = table;
count = 0;
for i = 1:nCells
    switch cell_types{i}
        case {'stellate','Stellate'}
            name = 'Stellate';
        case {'pyramidal','Pyramidal'}
            name = 'Pyramidal';
        case {'fast spiking','FastSpiking'}
            name = 'FastSpiking';
        otherwise
            name = cell_types{i};
    end
    if nGroups > 1
        for iGroup = 1:nGroups
            count = count + 1;
            fullname = [name '_' sprintf('group%d',iGroup)];
            data_table.(fullname) = data_all_mat(:,count);
        end
    else
        data_table.(name) = data_all_mat(:,i);
    end
end

end
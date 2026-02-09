function data_table = myLogSummaryTable(data_all,cell_types)
nCells = numel(cell_types);

% construct data array with NaNs
maxNumEl = max(cellfun(@numel,data_all(:)));
data_all_pad = cellfun(@(x){padarray(x(:),[maxNumEl-numel(x),0],NaN,'post')}, data_all(:)); % Pad each vector with NaN values to equate lengths
data_all_mat = cell2mat(data_all_pad'); 

% if size(data_all,2) == 1
%     data_all_mat = data_all_mat';
% end

% column names for cell types
varNames = cell(nCells,1);
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
    varNames{i} = name;
end

% Default: per column, 95% CI, log base = 10
out = geom_summary_array(data_all_mat, 'LogBase', 10, 'Conf', 0.95, 'Dim', 1, ...
                         'VarNames', varNames);

% View the tidy table:
data_table = out.Table;


% data_table = array2table([meanData; sd; SEM; N; medianData; Q(1,:); Q(2,:); dataMin; dataMax; dataRange],'VariableNames',varNames,'RowNames',{'Mean','SD','SEM','N','Median','Q25','Q75','Min','Max','Range'});

end
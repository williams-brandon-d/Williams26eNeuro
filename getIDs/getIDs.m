function IDs = getIDs(info,ps)
% find IDs for analysis
if strcmp(ps.locations,'all'); ps.locations = unique({info.location}); end
if strcmp(ps.cell_types,'all'); ps.cell_types = unique({info.cell_type}); end
if strcmp(ps.experiments,'all'); ps.experiments = unique({info.experiment}); end
if strcmp(ps.comments,'all'); ps.comments = unique({info.comments}); end
if strcmp(ps.protocols,'all'); ps.protocols = unique({info.protocol}); end
if strcmp(ps.cell_nums,'all')
    p.cell_num = 'all';
else
    p.cell_num = ps.cell_nums;
end

fieldElements = structfun(@numel,ps);
nCells = prod(fieldElements);
IDs_cell = cell(nCells,1);
i = 0;
for iLocation = 1:numel(ps.locations)
    p.location = ps.locations{iLocation};
    for iCellType = 1:numel(ps.cell_types)
        p.cell_type = ps.cell_types{iCellType};
        for iExperiment = 1:numel(ps.experiments)
            p.experiment = ps.experiments{iExperiment};
            for iComment = 1:numel(ps.comments)
                p.comments = ps.comments{iComment};
                for iProtocol = 1:numel(ps.protocols)
                    p.protocol = ps.protocols{iProtocol};
                    i = i + 1;
                    IDs_cell{i} = find_IDs(info,p); % append cell with IDs
                end
            end
        end
    end
end
IDs = vertcat(IDs_cell{:});

end
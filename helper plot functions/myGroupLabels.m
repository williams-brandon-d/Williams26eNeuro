function groupLabels = myGroupLabels(cell_types)

nCells = numel(cell_types);

groupLabels = cell(nCells,1);
for i = 1:nCells
    celltype = cell_types{i};
    switch celltype
        case 'stellate'
            name = 'SC';
        case 'pyramidal'
            name = 'PC';
        case 'fast spiking'
            name = 'FS';
    end
    groupLabels(i) = {sprintf('%s',name)};
end
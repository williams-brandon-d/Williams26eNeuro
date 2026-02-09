function [colors,alphas] = getColorAlpha(cell_types,nCommentgroups)

nCells = numel(cell_types);

colors = cell(nCells*nCommentgroups,1);
alphas = colors;

count = 0;
for i = 1:nCells
celltype = cell_types{i};
    for ii = 1:nCommentgroups
        count = count + 1;
        switch celltype
            case 'stellate'
                colors(count) = {[1 0 0]}; % red 
            case 'pyramidal'
                colors(count) = {[0.4660, 0.6740, 0.1880]}; % green
            case 'fast spiking'
                colors(count) = {[0 0 1]}; % blue
        end
        if ii == 1 
            alphas{count} = 0.7; % control
        elseif ii == 2
            alphas{count} = 0.3; % dnqx
        end

    end
end


end
function [color,alpha] = getCellColorAlpha(cellType,comment)


switch cellType
    case 'stellate'
        color = [1 0 0]; % red 
    case 'pyramidal'
        color = [0.4660, 0.6740, 0.1880]; % green
    case 'fast spiking'
        color = [0 0 1]; % blue
end


switch comment
    case {'DNQX after'}
        alpha = 0.3;
    otherwise
        alpha = 0.7;
end




end
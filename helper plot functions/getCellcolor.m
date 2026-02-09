function [color,alpha] = getCellcolor(cellType,dataSet,comment)

        cweight = 1;

switch comment
    case {'DNQX after', 'DNQX after 10 uM'}
        alpha = 1; % 0.7
    case {'DNQX before','DNQX before 10 uM'}
        alpha = 1; % 0.7
    otherwise
        alpha = 1;
end

switch cellType
    case 'stellate'
        switch dataSet
            case 'Thy1'
                color = cweight*[1 0 0];           
            case 'Camk2'
                color = cweight*[1 0 0];
            case 'PV Transgenic'
                color = cweight*[1 0 0];
%                 color = cweight*[0.6350 0.0780 0.1840];
            case 'PV Viral'
                color = cweight*[1 0 0];
        end
    case 'pyramidal'
        switch dataSet
            case 'Thy1'
                color = cweight*[0.4660, 0.6740, 0.1880];
            case 'Camk2'
                color = cweight*[0.4660, 0.6740, 0.1880];
            case 'PV Transgenic'
                color = cweight*[0.4660, 0.6740, 0.1880];
%                 color = cweight*[0.8500 0.2250 0.0980];
            case 'PV Viral'
                color = cweight*[0.4660, 0.6740, 0.1880];
        end
    case 'fast spiking'
        switch dataSet
            case 'Thy1'
                color = cweight*[0 0 1];
            case 'Camk2'
                color = cweight*[0 0 1];
            case 'PV Transgenic'
                color = cweight*[0 0 1];
%                 color = cweight*[0 0.4470 0.7410];
            case 'PV Viral'
                color = cweight*[0 0 1];
        end
    case 'slow interneuron'
        color = [0.4940 0.1840 0.5560];
end




end
function y = getYparams(dataType,protocol,axisType)

switch axisType

    case 'normal'

        switch dataType
            case 'stim'
                y.min = 0; y.max = 10; y.dy = 2;
                y.ticks = y.min:y.dy:y.max;
                y.tickLabels = compose('%d',y.ticks);
                y.labelstring = 'LED Input (mV)';   
                y.scale = 'linear';
            case 'phase' 
                switch protocol
                    case 'theta'
                        y.min = -pi; y.max = pi; y.dy = pi;
                        y.ticks = y.min:y.dy:y.max;
                        y.tickLabels = {'-π','0','π'};
                        y.labelstring = 'Theta Phase (rad)';
                    case 'pulse'
                        y.min = 0; y.max = 100; y.dy = 50;
                        y.ticks = y.min:y.dy:y.max;
                        y.tickLabels = compose('%d',y.ticks);
                        y.labelstring = 'Time (ms)';
                end
                y.scale = 'linear';
        
            case {'frequency','PSDfrequency','frequencyACF'}
                y.min = 0; y.max = 200; y.dy = 50;
                y.ticks = y.min:y.dy:y.max;
                y.tickLabels = compose('%d',y.ticks);
                y.labelstring = 'Frequency (Hz)';
                y.scale = 'linear';
        
            case {'power'}
                switch protocol
                    case 'theta'
                    y.min = -1; y.max = 5; % ymin = 10^-1 for NRSA
                    case 'pulse'
                    y.min = -1; y.max = 5; % ymin = 10^1 for NRSA
                end
    %             y.ticks = logspace(0,4,5);
                y.ticks = y.min:y.max;
                y.tickLabels = compose('%d',y.ticks);
                y.labelstring = 'Log Peak Gamma Power (pA^{2})';
                y.scale = 'linear';
            case {'totalFastPower','totalGammaPower'}
                switch protocol
                    case 'theta'
                    y.min = 3; y.max = 9; % ymin = 10^-1 for NRSA
                    case 'pulse'
                    y.min = 3; y.max = 9; % ymin = 10^1 for NRSA
                end
    %             y.ticks = logspace(0,4,5);
                y.ticks = y.min:y.max;
                y.tickLabels = compose('%d',y.ticks);
                switch dataType
                    case 'totalGammaPower'
                        y.labelstring = 'Log Total Gamma Power (pA^{2})';
                    case 'totalFastPower'
                        y.labelstring = 'Log Total Fast Power (pA^{2})';
                end
                y.scale = 'linear';
            case {'PSDGammaPower','PSDFastPower'}
                switch protocol
                    case 'theta'
                    y.min = -1; y.max = 5; % ymin = 10^-1 for NRSA
                    case 'pulse'
                    y.min = -1; y.max = 5; % ymin = 10^1 for NRSA
                end
    %             y.ticks = logspace(0,4,5);
                y.ticks = y.min:y.max;
                y.tickLabels = compose('%d',y.ticks);
                switch dataType
                    case 'PSDGammaPower'
                        y.labelstring = 'Log Total Power (pA^{2})';
                    case 'PSDFastPower'
                        y.labelstring = 'Log Total Power (pA^{2})';
                end
                y.scale = 'linear';
            case 'meanCorr'
                y.min = 0; y.max = 1; y.dy = 0.2;
                y.ticks = y.min:y.dy:y.max;
                y.tickLabels = compose('%.1f',y.ticks);
                y.labelstring = 'Correlation Coeff.';
                y.scale = 'linear';

            case 'meanCorrLag'
                y.min = -10; y.max = 10; y.dy = 2;
%                 y.min = -10; y.max = 10; y.dy = 2;
                y.ticks = y.min:y.dy:y.max;
                y.tickLabels = compose('%d',y.ticks);
                y.labelstring = 'Lag (ms)';
                y.scale = 'linear';

            case 'ISI_ratio'
                y.min = 0; y.max = 3; y.dy = 1;
                y.ticks = y.min:y.dy:y.max;
                y.tickLabels = compose('%d',y.ticks);
                y.labelstring = 'ISI Ratio 1/2';
                y.scale = 'linear';

            case 'sag_hyper'
                y.min = 0; y.max = 20; y.dy = 5;
                y.ticks = y.min:y.dy:y.max;
                y.tickLabels = compose('%d',y.ticks);
                y.labelstring = 'Sag hyperpolarizing (mV)';
                y.scale = 'linear';

            case 'sag_depol_thresh'
                y.min = 0; y.max = 10; y.dy = 2;
                y.ticks = y.min:y.dy:y.max;
                y.tickLabels = compose('%d',y.ticks);
                y.labelstring = 'Sag depolarizing (mV)';
                y.scale = 'linear';
            case 'dap'
                y.min = 0; y.max = 5; y.dy = 1;
                y.ticks = y.min:y.dy:y.max;
                y.tickLabels = compose('%d',y.ticks);
                y.labelstring = 'depolarizing After Potential (mV)';
                y.scale = 'linear';

            case 'latency'
                y.min = 0; y.max = 1000; y.dy = 200;
                y.ticks = y.min:y.dy:y.max;
                y.tickLabels = compose('%d',y.ticks);
                y.labelstring = 'Latency to first spike (ms)';
                y.scale = 'linear';

            case 'time_constant'
                y.min = 0; y.max = 60; y.dy = 10;
                y.ticks = y.min:y.dy:y.max;
                y.tickLabels = compose('%d',y.ticks);
                y.labelstring = 'Membrane Time Constant (ms)';
                y.scale = 'linear';
            case 'rmp'
                y.min = -90; y.max = -40; y.dy = 10;
                y.ticks = y.min:y.dy:y.max;
                y.tickLabels = compose('%d',y.ticks);
                y.labelstring = 'Resting Membrane Potential (mV)';   
                y.scale = 'linear';
        end

    case 'log'
        switch dataType
            case {'power','PSDpower'}
                switch protocol
                    case 'theta'
                    y.min = 10^-1; y.max = 10^5; % ymin = 10^-1 for NRSA
                    case 'pulse'
                    y.min = 10^-1; y.max = 10^5; % ymin = 10^1 for NRSA
                end
                y.ticks = logspace(-1,5,7);
    %             y.ticks = y.min:y.max;
                y.tickLabels = compose('%d',log10(y.ticks));
                y.labelstring = 'Log Peak Gamma Power (pA^{2})';
                y.scale = 'log';
        end

    case 'difference'

        switch dataType
            case 'phase' 
                switch protocol
                    case 'theta'
                        y.min = -2*pi; y.max = 2*pi; y.dy = 2*pi;
                        y.ticks = y.min:y.dy:y.max;
                        y.tickLabels = {'-2π','0','2π'};
                        y.labelstring = '\Delta Theta Phase (rad)';
                    case 'pulse'
                        y.min = 0; y.max = 100; y.dy = 50;
                        y.ticks = y.min:y.dy:y.max;
                        y.tickLabels = compose('%d',y.ticks);
                        y.labelstring = '\Delta Time (ms)';
                end
                y.scale = 'linear';
        
            case {'frequency','frequencyACF'}
                y.min = -100; y.max = 100; y.dy = 25;
                y.ticks = y.min:y.dy:y.max;
                y.tickLabels = compose('%d',y.ticks);
                y.labelstring = '\Delta Frequency (Hz)';
                y.scale = 'linear';
        
            case 'power'
                switch protocol
                    case 'theta'
                    y.min = -4; y.max = 2; % ymin = 10^-1 for NRSA
                    case 'pulse'
                    y.min = -4; y.max = 2; % ymin = 10^1 for NRSA
                end
    %             y.ticks = logspace(0,4,5);
                y.ticks = y.min:y.max;
                y.tickLabels = compose('%d',y.ticks);
                y.labelstring = '\Delta Log Peak Gamma Power (pA^{2})';
                y.scale = 'linear';
           case {'totalFastPower','totalGammaPower'}
                switch protocol
                    case 'theta'
                    y.min = -3; y.max = 2; % ymin = 10^-1 for NRSA
                    case 'pulse'
                    y.min = -3; y.max = 2; % ymin = 10^1 for NRSA
                end
    %             y.ticks = logspace(0,4,5);
                y.ticks = y.min:y.max;
                y.tickLabels = compose('%d',y.ticks);
                switch dataType
                    case 'totalGammaPower'
                        y.labelstring = 'Log Total Gamma Power (pA^{2})';
                    case 'totalFastPower'
                        y.labelstring = 'Log Total Fast Power (pA^{2})';
                end
                y.scale = 'linear';
        end

    case 'normalized'
        y.min = 0;
        y.max = 100;
        y.dy = 25;
        y.ticks = y.min:y.dy:y.max;
        y.tickLabels = compose('%d',y.ticks);
        y.scale = 'linear';   
        switch dataType
            case 'power'
                y.labelstring = 'Peak Gamma Power (norm)';
            case {'frequency','frequencyACF'}
                y.labelstring = 'Frequency (norm)';
            case 'phase'
                y.labelstring = 'Theta Phase (norm)';
            case 'totalGammaPower'
                y.labelstring = 'Total Gamma Power (norm)';
            case 'totalFastPower'
                y.labelstring = 'Total Fast Power (norm)';
        end

    case 'percentdiff'
        y.min = -100;
        y.max = 100;
        y.dy = 50;
        y.ticks = y.min:y.dy:y.max;
        y.tickLabels = compose('%d',y.ticks);
        y.scale = 'linear';   
        switch dataType
            case 'power'
                y.labelstring = '%\Delta Peak Gamma Power';
            case {'frequency','frequencyACF'}
                y.labelstring = '%\Delta Frequency';
            case 'phase'
                y.labelstring = '%\Delta Theta Phase';
            case 'totalGammaPower'
                y.labelstring = '%\Delta Total Gamma Power (pA^{2})';
            case 'totalFastPower'
                y.labelstring = '%\Delta Total Fast Power (pA^{2})';
        end
end
function tau = calculateTimeConstant(data_sweep,Fs,pulse_start_index,pulse_stop_index,plotnum,cellType)
    
    dt = 1/Fs;

    data_pulse = data_sweep(pulse_start_index:pulse_stop_index);
    data_after = data_sweep(pulse_stop_index:end);
    time_after = (0:length(data_after)-1)*dt;

    ss_pulse = median(data_pulse); % find steady state during pulse
    ss_after = median(data_after); % find steady state value after pulse
    
    ss_range = ss_pulse - ss_after;
    value_90 = ss_after + 0.9*ss_range;
    value_10 = ss_after + 0.1*ss_range;
    
%     fprintf('value_90 = %.3g\n',value_90)
    [~,max_after_index] = max(data_after);
    data_after_max = data_after(max_after_index:end);

    exp_start = find(data_after_max < value_90,1);  
    exp_end = find(data_after_max < value_10,1);
    
    data_exp = data_after_max(exp_start:exp_end); % data to fit exponential
    time_exp = (0:length(data_exp)-1)*dt;
    time_exp_ms = time_exp*1000;
    
    norm_exp = (data_exp - min(data_exp))./(max(data_exp)-min(data_exp));
    modelfun = @(b,x) exp(-x./b(1));
    
    % estimate time constant based on cell type hypothesis
    switch cellType
        case 'fast spiking'
            beta0 = 5; % ms
        otherwise % stellate, pyramidal
            beta0 = 10; % ms
    end

    mdl = fitnlm(time_exp_ms,norm_exp,modelfun,beta0);

%     modelfun = @(b,x) b(1).*exp(b(2).*x)+b(3);
%     beta0 = [ss_range -0.1 value_10];
%     mdl = fitnlm(time_exp_ms,data_exp,modelfun,beta0);
    
    c = table2array(mdl.Coefficients(:,1));
    
    rsquared = mdl.Rsquared.Adjusted;

%     exp_fit = c(1).*exp(c(2).*time_exp_ms)+c(3);
%     tau = -1/c(2);

    exp_fit = range(data_exp).*exp(-time_exp_ms./c(1)) + min(data_exp);
    tau = c(1);

%     fprintf('tau = %.3g ms, R^2 = %.3g\n',tau,rsquared)
%     fprintf('c1 = %.3g mV, tau = %.3g ms, c3 = %.3g mV\n',c(1),tau,c(3))
    
    if plotnum == 1
        figure;
        plot(time_after,data_after)

        figure;
        h = plot(time_exp_ms,data_exp,'-b',time_exp_ms,exp_fit,'--r');
        set(h,'Linewidth',2)
        xlabel('Time (ms)')
        ylabel('Voltage (mV)')
        title(sprintf('Fit: Tau = %.3g ms, R^2 = %.3g',tau,rsquared))
    end
    
end
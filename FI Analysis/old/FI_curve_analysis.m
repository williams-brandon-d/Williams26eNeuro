%% FI curve analysis

clear variables; close all; clc;
cd('C:\Users\brndn\OneDrive\Desktop\White Lab\Matlab Files\');
addpath(genpath('./')); % add all folders and subfolders in cd to path

savenum = 1; % 0 = don't save, 1 = save info

dataSets = {'Thy1'};

% compare mean ISI ratio from all current levels

params.protocols = {'FI'};
params.experiments = {'currentclamp'};

params.location = 'all';
params.cell_num = 'all';
 
params.cell_types = {'stellate','pyramidal'};
params.comments = {'','DNQX before'}; % combine comment data


current_sweeps = -200:25:525; % pA

printnum = 0;

% 0 = no plot, 1 = plot data
sag_plotnum = 0;
ISI_plotnum = 0;
max_depol_plotnum = 0;
max_latency_plotnum = 0; 
max_dap_plotnum = 0; % just max sweep
data_plotnum = 0;
tau_plotnum = 0;

dap_plotnum = 0; % each sweep
depol_sag_plotnum = 0; % each sweep

pulse_start_sec = 1; % pulse starts at 1 sec
pulse_stop_sec = 2;

sag_sweep_index = 1; % current = -200 pA

spike_threshold = 0; % mV
spike_dist_ms = 1; % ms

data_lpf_fc = 2000; % Hz cutoff freq for window data
data_lpf_order = 8;

plot_sweeps = [];
plotTitle = ''; % if empty, auto generates title

tic;
for iSet = 1:numel(dataSets)
dataSet = dataSets{iSet};

[info,info_fullname,data_path] = getInfo(dataSet);

IDs = getIDs(info,params);

% main

if (isempty(IDs)); disp('No Files Found.'); return; end

nIDs = length(IDs);

for i = 1:nIDs
    ID = IDs{i};
    ID_index = find_index(info,'ID',ID);
    filename = sprintf('%s.abf',ID);

    p = info(ID_index);
    p.dataSet = dataSet;

    % folder save path for results --> location,cell_type,cell_num,experiment,comment,protocol
    if isempty(p.comments)
        comment = 'No Comment';
    else
        comment = p.comments;
    end

    saveFolder = sprintf('%sresults\\%s\\%s\\%s\\%s\\%s\\%s\\%s',data_path,p.location,p.cell_type,p.cell_num,p.experiment,p.protocol,comment,ID);
    
    fullDataName = fullfile(data_path,filename);
    
    % analyze file data 
    [data,si,file_info] = abfload(fullDataName,'start',0,'stop','e');
    data_units = char(file_info.recChUnits(1));
    fprintf('Analyzing %s,File %d/%d\n',filename,i,nIDs)

    dt = si*(1e-6); % sampling interval (seconds)
    Fs = 1/dt; % sampling frequency (Hz)
    
    data = squeeze(data);
    nSweeps = size(data,2);
        
    [b_lpf,a_lpf] = butter(data_lpf_order, data_lpf_fc/(Fs/2), 'low');
    data = filtfilt(b_lpf, a_lpf, data); 
    
    [sag_hyper,sag_depol] = calculateSAG(data,Fs,sag_plotnum);

    pulse_start_index = Fs*pulse_start_sec + 1;
    pulse_stop_index = Fs*pulse_stop_sec;
    
    spike_dist_sec = spike_dist_ms/1000;
    spike_point_distance = Fs*spike_dist_sec;

    sag_depol_thresh = cell(nSweeps,1);
    latency = cell(nSweeps,1);
    average_dap = cell(nSweeps,1);
    ISI_ratio = cell(nSweeps,1);

    dap_count = 0;
    perithreshold_sweep = [];

    for iSweep = 1:nSweeps
        data_sweep = data(:, iSweep);
        
        if (max(data_sweep) > spike_threshold) % prominence and width instead of threshold?
            [spike_pks, spike_indices] = findpeaks(data_sweep,'MinPeakHeight',spike_threshold,'MinPeakDistance',spike_point_distance);
            nSpikes = length(spike_indices);
            
            if current_sweeps(iSweep) < 300 && dap_count < 4
                [average_dap{iSweep},~] = calculateDAP(data_sweep,spike_indices,Fs,pulse_start_index,pulse_stop_index,dap_plotnum);
                if ~isempty(average_dap{iSweep}); dap_count = dap_count + 1; end
            end
            
            if (spike_indices(1) < pulse_stop_index) && (spike_indices(1) > pulse_start_index)
                latency{iSweep} = (spike_indices(1)-pulse_start_index)*dt*1000;
            end
            
            if nSpikes >= 3
                ISI_ratio{iSweep} = calculateISI(data_sweep,spike_indices,Fs,ISI_plotnum);
            end 
            
        elseif current_sweeps(iSweep) > 0
            sag_depol_thresh{iSweep} = calculateDepolSAG(data_sweep,Fs,pulse_start_index,pulse_stop_index,depol_sag_plotnum);
            
            if max(data_sweep(pulse_stop_index:end)) < spike_threshold % check for spikes after pulse
                perithreshold_sweep = iSweep;
            end
        end
    end
    
   if ~isempty(perithreshold_sweep)
      info(ID_index).time_constant = calculateTimeConstant(data(:,perithreshold_sweep),Fs,pulse_start_index,pulse_stop_index,tau_plotnum);
   end
    
    depol_not_empty = find(~cellfun(@isempty,sag_depol_thresh));
    sag_depol_mat = cell2mat(sag_depol_thresh);
    [sag_depol_threshold,max_depol_mat_index] = max(sag_depol_mat);
    max_depol_index = depol_not_empty(max_depol_mat_index);
    if max_depol_plotnum == 1 && ~isempty(max_depol_index)
        data_sweep = data(:,max_depol_index);
        calculateDepolSAG(data_sweep,Fs,pulse_start_index,pulse_stop_index,max_depol_plotnum);
    end
    
    latency_not_empty = find(~cellfun(@isempty,latency));
    latency_mat = cell2mat(latency);
    [max_latency,max_latency_mat_index] = max(latency_mat);
    max_latency_index = latency_not_empty(max_latency_mat_index);
    if max_latency_plotnum == 1
       data_sweep = data(:,max_latency_index);
       [~, spike_indices] = findpeaks(data_sweep,'MinPeakHeight',spike_threshold,'MinPeakDistance',spike_point_distance);
       plotMaxLatency(data_sweep,spike_indices,Fs)
    end
    
    dap_not_empty = find(~cellfun(@isempty,average_dap));
    dap_mat = cell2mat(average_dap);
    [max_dap,max_dap_mat_index] = max(dap_mat);
    if max_dap == 0
        max_dap_index = max_latency_index;
    else
        max_dap_index = dap_not_empty(max_dap_mat_index);
    end
    if max_dap_plotnum == 1
        data_sweep = data(:,max_dap_index);
        [~, spike_indices] = findpeaks(data_sweep,'MinPeakHeight',spike_threshold,'MinPeakDistance',spike_point_distance);
        [max_dap_cal,dap_cal] = calculateDAP(data_sweep,spike_indices,Fs,pulse_start_index,pulse_stop_index,max_dap_plotnum);
    end
    
    if data_plotnum == 1 
        if isempty( plot_sweeps )
            plot_sweeps = 1:ISI_sweep_index(1);
        elseif strcmp(plot_sweeps,'all')
            plot_sweeps = 1:nSweeps;            
        elseif max(plot_sweeps) > nSweeps
            if mod(nSweeps,2) == 1
                plot_sweeps = 1:2:nSweeps;
            else
                plot_sweeps = 1:2:(nSweeps-1);
            end
        end
        
        plotData(data,plot_sweeps,dt,plotTitle,info,ID_index)
    end 
    
    ISI_array = cell2mat(ISI_ratio);
    
    if printnum == 1
        fprintf('Hyper Sag = %.2f mV\n',sag_hyper)
        fprintf('Depol Sag = %.2f mV\n',sag_depol)
        fprintf('Threshold Depol Sag = %.2f mV\n',sag_depol_threshold)
        fprintf('First ISI 1/2 = %.2f\n',ISI_array(1))
        fprintf('Second ISI 1/2 = %.2f\n',ISI_array(2))
        fprintf('Latency to First Spike = %.2f ms\n',max_latency)
        fprintf('dAP = %.2f mV\n',max_dap)
    end
    
    info(ID_index).sag_depol = sag_depol;
    if sag_depol_threshold < 0
        sag_depol_threshold = 0;
    end
    info(ID_index).sag_depol_thresh = sag_depol_threshold;
    info(ID_index).sag_hyper = sag_hyper;
    info(ID_index).dap = max_dap;
    info(ID_index).latency = max_latency;
    
%     if ISI_array(1) > 2
%         info(ID_index).ISI_ratio = ISI_array(2);
%     else
%         info(ID_index).ISI_ratio = ISI_array(1);
%     end
    
    info(ID_index).ISI_ratio = median(ISI_array(1:3));
%     fprintf('Median ISI 1/2 = %.2f\n',info(ID_index).ISI_ratio)
%         fprintf('dAP = %.2f mV\n',max_dap)
end

if savenum == 1
    save(fullfile(save_path,save_name),'info')
end
toc;

end
%% functions

function [average_dap,dap] = calculateDAP(data_sweep,spike_indices,Fs,pulse_start_index,pulse_stop_index,plotnum)
    ISI_min = 20; % ms
    stop_time = 6; % ms
    stop_time_after_spike = stop_time/1000; % sec
    xlimits = [500 2500];

    nSpikes = length(spike_indices);
    nSamples = length(data_sweep);
    dt = 1/Fs;
    time = (1:nSamples)*dt*1000;
%     dap = [];
    dap = cell(nSpikes,1);
    plot_flag = 0;
    
    for iSpike = 1:nSpikes
        start_index = spike_indices(iSpike);
        
        if (start_index > pulse_stop_index) || (start_index < pulse_start_index)
            continue; % ignore spikes outside pulse
        end
        
        if iSpike < nSpikes
            ISI = (spike_indices(iSpike+1) - spike_indices(iSpike))*dt*1000;
            if ISI < ISI_min
                continue; % ignore small ISIs
            end
        end
        
        stop_index_after_spike = start_index + Fs*stop_time_after_spike;
        window_after_spike = data_sweep(start_index:stop_index_after_spike);
        
%         if iSpike == 1
%             dv = diff(window_after_spike);
%             dv2 = diff(dv);
%             figure;
%             plot(dv)
%             figure;
%             plot(dv2)
%         end
        
        [fAHP,fAHP_index] = min(window_after_spike);
        
        dap_window = window_after_spike(fAHP_index:end);
        [max_after_min,MaM_index] = max(dap_window);
        
        dap_iSpike = max_after_min - fAHP;
        
        if dap_iSpike > 30
           dap_iSpike = []; % ignore spikes ?
        end
        
%         dap = [dap dap_iSpike];
        dap{iSpike} = dap_iSpike;
        
        plot_flag = plot_flag + 1;
        
        if plotnum == 1 && plot_flag == 1
            figure;
            plot(time,data_sweep,'-k')
            hold on
            plot(time(start_index+fAHP_index-1),fAHP,'xb')
            plot( [time(1) time(end)],[fAHP fAHP],'--b' )
            plot(time(start_index+fAHP_index+MaM_index-2),max_after_min,'ob')
            plot( [time(1) time(end)],[max_after_min max_after_min],'--b' )
            plot( [time(start_index+fAHP_index+MaM_index-2) time(start_index+fAHP_index+MaM_index-2)],[fAHP max_after_min],'--k' )
            title( sprintf( 'dAP = %.1f mV', dap_iSpike ) )
            if ~isempty(xlimits); xlim(xlimits); end
            xlabel('Time (ms)')
            ylabel('Membrane Voltage (mV)')
        end
        
    end
    dap_mat = cell2mat(dap);
    average_dap = mean(dap_mat);
end




function [sag_hyper,sag_depol] = calculateSAG(data,Fs,plotnum)

    pulse_start_sec = 1; % pulse starts at 1 sec
    pulse_stop_sec = 2; % pulse ends at 2 sec
    sag_sweep_index = 1; % 1st current sweep, -200 pA
    
    pulse_start_index = Fs*pulse_start_sec + 1;
    pulse_stop_index = Fs*pulse_stop_sec;

    data_hyper_sag = data( pulse_start_index:pulse_stop_index, sag_sweep_index );
    pulse_ss = median( data_hyper_sag ); % steady state of sweep during pulse
    [pulse_min,min_index] = min( data_hyper_sag ); % peak hyper sag
    sag_hyper = pulse_ss - pulse_min;
    
    depol_start_index = pulse_stop_index + 1;
    data_depol_sag = data( depol_start_index:end, sag_sweep_index );
    after_pulse_ss = median( data_depol_sag ); % steady state after pulse
    [after_pulse_max,~] = max( data_depol_sag ); % peak depol sag
    sag_depol = after_pulse_max - after_pulse_ss;
    
    if sag_depol > 50; sag_depol = []; end % ignore spike
    
    if plotnum == 1
        nSamples = size(data,1);
        dt = 1/Fs;
        time = (1:nSamples)*dt*1000;
        time_hyper_peak = time( pulse_start_index + min_index );
%         time_depol_peak = time( depol_start_index + max_index );

        figure;
        plot(time, data(:,sag_sweep_index),'-k','Linewidth', 2)
        hold on
        plot([time(1) time(end)],[pulse_ss pulse_ss],'--r')
        plot([time(1) time(end)],[pulse_min pulse_min],'--r')
        plot([time_hyper_peak time_hyper_peak],[pulse_min pulse_ss],'--r')
        plot(time_hyper_peak,pulse_min,'or')
%         if ~isempty( sag_depol )
%             plot([time(1) time(end)],[after_pulse_ss after_pulse_ss],'--b')
%             plot([time(1) time(end)],[after_pulse_max after_pulse_max],'--b')
%             plot([time_depol_peak time_depol_peak],[after_pulse_max after_pulse_ss],'--b')
%             plot(time_depol_peak,after_pulse_max,'ob')
%         end
        xlabel('Time (ms)')
        ylabel('Membrane Voltage (mV)')
        title( sprintf('Hyper. SAG = %.1f mV',sag_hyper) )
%         title( sprintf('Hyper. SAG = %.1f mV, Depol. SAG = %.1f mV',sag_hyper,sag_depol) )
    end
end




function ISI_ratio = calculateISI(data_sweep,spike_indices,dt,plotnum)
        
    ISI_1 = (spike_indices(2) - spike_indices(1))*dt*1000;
    ISI_2 = (spike_indices(3) - spike_indices(2))*dt*1000;
    ISI_ratio = ISI_1/ISI_2;
    
    if plotnum == 1
        xlim_space = 5; % ms
        spike_line_min = -100; % mV
        spike_line_max = 100; % mV
        ISI_line_voltage = 0;
        ISI_color = {'b','g'};
        
        nSamples = size(data_sweep,1);
        time = (1:nSamples)*dt*1000;
        
        xlim_start = time(spike_indices(1)) - xlim_space;
        xlim_stop = time(spike_indices(3)) + xlim_space;
        xlimits = [xlim_start xlim_stop];

        figure;
        plot(time,data_sweep,'-k') % plot data
        hold on
        plot(time(spike_indices),data_sweep(spike_indices),'or') % plot spike peaks
        for iSpike = 1:3
            spike_index = spike_indices(iSpike);
            time_spike_index = time(spike_index);
            plot([time_spike_index time_spike_index],[spike_line_min spike_line_max],'--r') % plot spike time lines
        end
        for iISI = 1:2
            time_spike1 = time(spike_indices(iISI));
            time_spike2 = time(spike_indices(iISI+1));
            plot([time_spike1 time_spike2],[ISI_line_voltage ISI_line_voltage],'Linestyle','--','Color',ISI_color{iISI})
        end
        xlim(xlimits)
        xlabel('Time (ms)')
        ylabel('Membrane Voltage (mV)')
        title( sprintf( 'ISI 1/2 = %.2f', ISI_ratio) )
    end
end





function sag_depol = calculateDepolSAG(data_sweep,Fs,pulse_start_index,pulse_stop_index,plotnum)

    depol_stop_sec = 1.2;
    depol_stop_index = Fs*depol_stop_sec;
    
    data_pulse = data_sweep( pulse_start_index:pulse_stop_index );
    pulse_ss = median( data_pulse ); % steady state after pulse
    
    data_depol = data_sweep(pulse_start_index:depol_stop_index);
    [depol_max,max_index] = max( data_depol ); % peak depol sag
    sag_depol = depol_max - pulse_ss;
    
    if plotnum == 1
        nSamples = size(data_sweep,1);
        dt = 1/Fs;
        time = (1:nSamples)*dt*1000;
        time_depol_peak = time( pulse_start_index + max_index );

        figure;
        plot(time, data_sweep,'-k','Linewidth', 2)
        hold on
        plot([time(1) time(end)],[pulse_ss pulse_ss],'--b')
        plot([time(1) time(end)],[depol_max depol_max],'--b')
        plot([time_depol_peak time_depol_peak],[depol_max pulse_ss],'--b')
        plot(time_depol_peak,depol_max,'ob')
        xlabel('Time (ms)')
        ylabel('Membrane Voltage (mV)')
        title( sprintf('Depol. SAG = %.1f mV',sag_depol) )
    end
end





function plotMaxLatency(data_sweep,spike_indices,Fs)
    xlimits = [500 2500]; % Time (ms) for viewing only not analysis
    min_line_voltage = -100;
    max_line_voltage = 100;
    connect_line_voltage = 0;
    pulse_start_time = 1;

    dt = 1/Fs;
    pulse_start_index = Fs*pulse_start_time;
    
    nSamples = size(data_sweep,1);
    time = (0:nSamples-1)*dt*1000; % time in msec
    time = time';
    
    first_spike_time_ms = time(spike_indices(1));
    pulse_start_time_ms = time(pulse_start_index);
    max_latency = first_spike_time_ms - pulse_start_time_ms;

    figure;
    plot(time,data_sweep,'-k')
    hold on
    if ~isempty(xlimits); xlim(xlimits); end
    plot([pulse_start_time_ms pulse_start_time_ms],[min_line_voltage max_line_voltage],'--r')
    plot([first_spike_time_ms first_spike_time_ms],[min_line_voltage max_line_voltage],'--r')
    plot([pulse_start_time_ms first_spike_time_ms],[connect_line_voltage connect_line_voltage],'--b')
    xlabel('Time (ms)')
    ylabel('Membrane Voltage (mV)')
    title( sprintf('Max Latency = %.2f ms',max_latency) )

end





function plotData(data,plot_sweeps,dt,plotTitle,info,ID_index)

    xlimits = [500 2500]; % Time (ms) for viewing only not analysis
    ylimits = [];
    plotTimeOffset = 1000;
    plotLinewidth = 1;
    tickfontsize = 15;
    labelfontsize = 20;
    titlefontsize = 25;

    nSamples = size(data,1);
    time = (0:nSamples-1)*dt*1000 - plotTimeOffset; % time in msec
    time = time';

    figure;
    plot(time,data(:,plot_sweeps),'-k','Linewidth',plotLinewidth)
    if ~isempty(xlimits); xlim(xlimits-plotTimeOffset); end
    if ~isempty(ylimits); ylim(ylimits); end
    xlabel('Time (ms)','FontSize',labelfontsize,'FontWeight','bold')
    ylabel('Membrane Voltage (mV)','FontSize',labelfontsize,'FontWeight','bold')
    ax = gca;
    ax.YAxis.FontSize = tickfontsize;
    ax.XAxis.FontSize = tickfontsize;
    ax.YAxis.FontWeight = 'bold';
    ax.XAxis.FontWeight = 'bold';
    box off
%         set(gca,'Visible','off') % turn everything off besides trace
    if isempty(plotTitle)
       han = title(strjoin({info(ID_index).location info(ID_index).cell_type info(ID_index).cell_num }),'FontSize',titlefontsize,'FontWeight','bold');
    else
       han = title(plotTitle,'FontSize',titlefontsize,'FontWeight','bold');
    end
    set(han,'Visible','on') % make title visible
end

function tau = calculateTimeConstant(data_sweep,Fs,pulse_start_index,pulse_stop_index,plotnum)
    
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
    beta0 = 10; % ms
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
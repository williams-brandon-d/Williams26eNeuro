%% FI curve analysis
clear variables
close all
clc

info_path = 'C:\Users\brndn\OneDrive\Desktop\White Lab\thy1_chr2\Thy1-ChR2\MATLAB\thy1chr2\';
info_name = 'thy1_chr2_121920.mat';

data_path = 'C:\Users\brndn\Downloads\Thy1-ChR2\Raw Data\mEC\';

save_path = 'C:\Users\brndn\OneDrive\Desktop\White Lab\My Papers\Paper 1\Figures\Supp Fig 1\';

% Thy1
% Ex) dorsal stellate 12 sweeps = 1:2:19;
% Ex) ventral pyramidal 5 sweeps = 1:2:13;
% Ex) dorsal fast spiking 3 sweeps = 1:2:15;

p.location = 'dorsal';
p.cell_type = 'fast spiking';
p.cell_num = '3';
p.comments = '';

p.protocol = 'FI';
p.experiment = 'currentclamp';

save_num = 1; % 0 = don't save info, 1 = save info
data_plotnum = 1; % 0 = no plot, 1 = plot data

xlimits = [750 2250]; % Time (ms) for viewing only not analysis
ylimits = [-100 50];

switch p.cell_type
    case 'stellate'
        plot_sweeps = 1:2:19;
        plotTitle = 'Stellate'; % if empty, auto generates title
    case 'pyramidal'
        plot_sweeps = 1:2:13;
        plotTitle = 'Pyramidal'; % if empty, auto generates title
    case 'fast spiking'
        plot_sweeps = 1:2:15;
        plotTitle = 'Fast Spiking'; % if empty, auto generates title
end

save_name = sprintf('%s example',plotTitle);

plotTimeOffset = 1000;
plotLinewidth = 1;
tickfontsize = 15;
labelfontsize = 20;
titlefontsize = 25;



load(fullfile(info_path,info_name),'info')

IDs = find_IDs(info,p); % find IDs based on parameters in "p"

if (isempty(IDs)); disp('No Files Found.'); return; end

for i = 1:length(IDs)
    ID = IDs{i};
    ID_index = find_index(info,'ID',ID);
    filename = sprintf('%s.abf',ID);
    fullname = fullfile(data_path,filename);
    [data,si,file_info] = abfload(fullname,'start',0,'stop','e');
    
    data_units = char(file_info.recChUnits(1));
    
    dt = si*(1e-6); % sampling interval (seconds)
    Fs = 1/dt; % sampling frequency (Hz)
    
    data = squeeze(data);
    [nSamples,nSweeps] = size(data);
    
    time = (0:nSamples-1)*si*(1e-3) - plotTimeOffset; % time in msec
    time = time';
    
    if max(plot_sweeps) > nSweeps
        if mod(nSweeps,2) == 1
            plot_sweeps = 1:2:nSweeps;
        else
            plot_sweeps = 1:2:(nSweeps-1);
        end
    end
    
    data_sweeps = data(:,plot_sweeps);
  
    data_lpf_fc = 2000; % Hz cutoff freq for window data
    data_lpf_order = 8;
    
    [b_lpf,a_lpf] = butter(data_lpf_order,data_lpf_fc/(Fs/2),'low');
    lpf_data = filtfilt(b_lpf,a_lpf,data_sweeps);
    
    if data_plotnum == 1
        [colors,alphas] = getColorAlpha({p.cell_type},1);
        figure;
        plot(time,lpf_data,'-','Linewidth',plotLinewidth,'Color',colors{1})
        xlim(xlimits-plotTimeOffset)
        ylim(ylimits)
        xlabel('Time (ms)','FontSize',labelfontsize,'FontWeight','bold')
        ylabel(['Membrane Voltage (' data_units ')'],'FontSize',labelfontsize,'FontWeight','bold')
        ax = gca;
        ax.YAxis.FontSize = tickfontsize;
        ax.XAxis.FontSize = tickfontsize;
        ax.YAxis.FontWeight = 'bold';
        ax.XAxis.FontWeight = 'bold';
        box off
%         set(gca,'Visible','off')
        if isempty(plotTitle)
            han = title(strjoin({info(ID_index).location info(ID_index).cell_type info(ID_index).cell_num }),'FontSize',titlefontsize,'FontWeight','bold');
        else
           han = title(plotTitle,'FontSize',titlefontsize,'FontWeight','bold');
        end
        set(han,'Visible','on')
        if save_num == 1
            saveas(gcf,fullfile(save_path,save_name),'svg') 
        end
    end  
end
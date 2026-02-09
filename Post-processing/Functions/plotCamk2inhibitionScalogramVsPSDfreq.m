function plotCamk2inhibitionScalogramVsPSDfreq()

clear variables; close all; clc;
cd('C:\Users\brndn\OneDrive\Desktop\White Lab\Matlab Files\');

savenum = 1; % 0 = don't save, 1 = save info

% plot inhibition vs LFP frequency

% for duplicates - take largest value for inhibition?
% should be no duplicates anymore


dataSets = {'Camk2'};
cell_types = {'stellate','pyramidal','fast spiking'};

% find ID parameters for analysis
params.locations = 'all';
params.experiments = {'inhibition'};
params.cell_nums = 'all';
params.comments = {'','Gabazine before'}; % 'GBZ before' 'DNQX_GABAzine_before', 'Gabazine before', 'DNQX before'
params.protocols = {'theta_2chan'};

% if strcmp(dataSets,'all'); dataSets = {'Camk2','Thy1','PV Transgenic','PV Viral'}; end

% load data for figure

nCellTypes = numel(cell_types);

for iSet = 1:numel(dataSets)
dataSet = dataSets{iSet};

[info,info_fullname,data_path] = getInfo(dataSet);

xCell = cell(nCellTypes,1);
yCell = cell(nCellTypes,1);

for iCell = 1:nCellTypes
params.cell_types = cell_types(iCell);
IDs = getIDs(info,params);

validIDs = removeBadLFP(IDs,info,dataSet); % remove bad LFP recordings (CamK2 only)

if (isempty(IDs)); fprintf('No %s Files Found.',params.cell_types); continue; end

nIDs = numel(IDs);

xArray = zeros(nIDs,1);
yArray = zeros(nIDs,1);

for iID = 1:nIDs

    ID = IDs{iID};

    % ignore bad IDs
    if ~any(strcmp(ID,validIDs))
        xArray(iID) = NaN; % ignore if either value is empty
        yArray(iID) = NaN;
        continue;
    end

    ID_index = find_index(info,'ID',ID);
    filename = sprintf('%s.abf',ID);
    fprintf('Analyzing %s,File %d/%d\n',filename,iID,nIDs)
    p = info(ID_index);

    if isempty(p.comments)
        comment = 'No Comment';
    else
        comment = p.comments;
    end

    commentFolder = sprintf('%sresults\\%s\\%s\\%s\\%s\\%s\\%s\\%s',data_path,p.location,p.cell_type,p.cell_num,p.experiment,p.protocol,comment);

    dataFolder = getIDFolder(commentFolder,ID);

    dataFilename = [dataFolder filesep 'data.mat'];

    load(dataFilename,'file');

    % gather relevant data 
    if ~isempty(file.cell.CWTmaxValues) && ~isempty(file.lfp.CWTmaxValues)
        xArray(iID) = file.cell.CWTmaxValues(2); % cell peak scalogram freq
        yArray(iID) = file.cell.psd.max_psd_gamma_freq; % cell peak psd freq
%         yArray(iID) = file.cell.psdAll.max_psd_gamma_freq_all; % cell peak psd freq
    else
        xArray(iID) = NaN; % ignore if either value is empty
        yArray(iID) = NaN;
    end

    % add same exclusion criteria as boxplots
    removeFlag = removeLowPowerData(file,dataSet);
    if removeFlag; xArray(iID) = NaN; yArray(iID) = NaN; end

    clearvars file;
end

xCell{iCell} = xArray;
yCell{iCell} = yArray;

end

end

% vectorize cell data
x = rmmissing(vertcat(xCell{:}));
y = rmmissing(vertcat(yCell{:}));

% linear fit data
[yfit,Rsq] = linearfit(x,y);

% plot figure
fig = figure;
hold on;
for iCell = 1:nCellTypes
    cellType = cell_types{iCell};

    switch cellType
        case 'stellate'
            color = [1 0 0];
        case 'pyramidal'
            color = [0 1 0];
        case 'fast spiking'
            color = [0 0 1];
    end

    plot(xCell{iCell},yCell{iCell},'.','Color',color,'DisplayName',cellType)

end

han1 = plot(x,yfit,'-k','Linewidth',2,'DisplayName',sprintf('   R^{2} = %.3f',Rsq)); 

fontsize = 12;
fontweight = 'bold';
ax = gca;
ax.YAxis.FontSize = fontsize;
ax.XAxis.FontSize = fontsize;
ax.YAxis.FontWeight = fontweight;
ax.XAxis.FontWeight = fontweight;

axis square;

limits = [50 150]; % Hz
xlim(limits);
ylim(limits);

xlabel('CaMK2 Inhibition Scalogram Freq (Hz)','FontSize',12,'FontWeight','bold');
ylabel('CaMK2 Inhibition PSD Freq (Hz)','FontSize',12,'FontWeight','bold');

legend(han1,'location','northeastoutside','FontSize',12,'FontWeight','bold')
legend('boxoff');

if savenum
    saveFolder = 'C:\Users\brndn\OneDrive\Desktop\White Lab\Carmen Grants';
    saveFilename = [saveFolder filesep 'camk2 inhibition scalogram vs psd freq.svg'];
    print(fig,'-vector','-dsvg',saveFilename);
end

close all;

    function [yfit,Rsq] = linearfit(x,y)
        % x and y are vectors
        x = reshape(x,[],1); % column vector
        y = reshape(y,[],1); % column vector
        % linear fit for all data
        X = [ones(length(x),1) x];
        b = X \ y;
        yfit = X*b; % linear fit
        Rsq = 1 - sum((y - yfit).^2)/sum((y - mean(y)).^2);

    end

end
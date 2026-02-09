function plotLightvsCell()

clear variables; close all; clc;
cd('C:\Users\brndn\OneDrive\Desktop\White Lab\Matlab Files\');

savenum = 1; % 0 = don't save, 1 = save info


% plot light power vs gamma power/freq, average firing rate

% test nonlinear monotonic correlation - spearman??
% [RHO,PVAL] = corr(a',b','Type','Spearman');

% convert light input mV to light intensity (mW/mm^2)

dataSets = {'Thy1','PV Transgenic'};
cell_types = {'stellate','pyramidal','fast spiking'};
% cell_types = {'stellate','pyramidal'};

dataTypes = {'power','freq','rate'};

% find ID parameters for analysis
params.locations = 'all';
params.cell_nums = 'all';
params.comments = {'','DNQX before','DNQX before 10 uM','GBZ before','DNQX_GABAzine_before', 'Gabazine before'}; % 'GBZ before' 'DNQX_GABAzine_before', 'Gabazine before', 'DNQX before'
params.protocols = {'theta','theta_2chan'};

if strcmp(dataSets,'all'); dataSets = {'Camk2','Thy1','PV Transgenic','PV Viral'}; end


% load data for figure

nCellTypes = numel(cell_types);
nDataTypes = numel(dataTypes);

for iSet = 1:numel(dataSets)
dataSet = dataSets{iSet};

[info,info_fullname,data_path] = getInfo(dataSet);

switch dataSet
    case 'Thy1'
         saveFolder = 'C:\Users\brndn\Downloads\Thy1-ChR2\Raw Data\mEC\results\Summary';
%          experiments = {'inhibition'};
    case 'PV Transgenic'
         saveFolder = 'C:\Users\brndn\Downloads\PV-ChR2 Transgenic\Summary';
    case 'PV Viral'
         saveFolder = 'C:\Users\brndn\Downloads\PV-ChR2\Summary';
    case 'Camk2'
         saveFolder = 'C:\Users\brndn\Downloads\CaMK2-ChR2\Summary';
%          experiments = {'excitation'};
end

for iType = 1:nDataTypes
dataType = dataTypes{iType};

switch dataType
    case 'power'
        experiment = 'inhibition';
    case 'freq'
        experiment = 'inhibition';
    case 'rate'
        experiment = 'currentclamp';
end

params.experiments = {experiment};

xCell = cell(nCellTypes,1);
yCell = cell(nCellTypes,1);

for iCell = 1:nCellTypes
params.cell_types = {cell_types{iCell}};
IDs = getIDs(info,params);

IDs = removeIDs(IDs,info,dataSet); % skip bad recordings - get params for cells to skip

% validIDs = removeBadLFP(IDs,info,dataSet); % remove bad LFP recordings (CamK2 only)

if (isempty(IDs)); fprintf('No %s Files Found.',cell_types{iCell}); continue; end

nIDs = numel(IDs);

inputArray = zeros(nIDs,1);
xArray = zeros(nIDs,1);
yArray = zeros(nIDs,1);

for iID = 1:nIDs

    ID = IDs{iID};

%     % ignore bad IDs
%     if ~any(strcmp(ID,validIDs))
%         xArray(iID) = NaN; % ignore if either value is empty
%         yArray(iID) = NaN;
%         continue;
%     end

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
    if ~isempty(file.led_input)
        focalArea = 0.237; % mm^2 40x objective 
        % convert led input to light intensity
        power = LEDvoltageToPower470nm(file.led_input); % mW
        intensity = power / focalArea; % mW / mm^2
        inputArray(iID) = file.led_input;
        xArray(iID) = intensity;
        switch dataType 
            case 'power'
                yArray(iID) = log10(file.cell.CWTmaxValues(3)); % log lfp peak scalogram power
            case 'freq'
                yArray(iID) = file.cell.CWTmaxValues(2); % lfp peak scalogram freq
            case 'rate'
                % calculate average firing rate
%                 file.FRs = FRs;
                nSpikes = numel([file.theta_spike_phases{:}]); % nCycles x cycle phase array
                nCycles = numel(file.cycles);
                if ~isempty(nSpikes)
                    yArray(iID) = nSpikes / nCycles; % average spikes per theta cycle 
                else
                    yArray(iID) = 0; % change to NaN to ignore non-firing cells
                end
        end
    else
        xArray(iID) = NaN; % ignore if either value is empty
        yArray(iID) = NaN;
    end

    clearvars file;
end

xCell{iCell} = xArray;
yCell{iCell} = yArray;

end

% vectorize cell data
x = rmmissing(vertcat(xCell{:}));
y = rmmissing(vertcat(yCell{:}));

% save all data
saveFilename = [saveFolder filesep sprintf('Light vs %s stats.xlsx',dataType)];
stats = myMultipleIndependentGroupStats({x,y},{'Light Intensity',dataType},'','linear',saveFilename,1);

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


xlimits = [0 25]; % mW / mm^2
xlim(xlimits);
xlabel('Light Intensity (mW/mm^{2})','FontSize',12,'FontWeight','bold');

switch dataType
    case 'power'
        ylimits = [0 5]; % log power
        ystring = sprintf('%s Peak Gamma Power (Log pA^{2})',dataSet);
    case 'freq'
        ylimits = [50 200]; % Hz
        ystring = sprintf('%s Peak Frequency (Hz)',dataSet);
    case 'rate'
        ylimits = [0 5]; % Hz
        ystring = sprintf('%s Firing Rate (Spikes / Theta Cycle)',dataSet);      
end

ylim(ylimits);
ylabel(ystring,'FontSize',12,'FontWeight','bold');

legend(han1,'location','northeastoutside','FontSize',12,'FontWeight','bold')
legend('boxoff');

if savenum
%     saveFolder = 'C:\Users\brndn\OneDrive\Desktop\White Lab\Carmen Grants';
    saveFilename = [saveFolder filesep sprintf('%s light power vs %s .svg',dataSet,dataType)];
    print(fig,'-vector','-dsvg',saveFilename);
end

end

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

    function powers = LEDvoltageToPower470nm(voltages)
        % measured LED input voltages and light powers
        measuredVoltages =    [0.500 0.600 0.700 0.800 0.900 1.000 1.500 2.000 5.000 10.000]; % mV
        measuredPowers = [0.000 0.040 0.078 0.110 0.145 0.180 0.350 0.495 1.425 2.8450]; % mW

        % linear interpolation
        powers = interp1(measuredVoltages,measuredPowers,voltages,"linear",'extrap');

    end

end
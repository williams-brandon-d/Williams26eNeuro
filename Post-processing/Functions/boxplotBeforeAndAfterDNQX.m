function boxplotBeforeAndAfterDNQX(dataSets, comment_group1,comment_group2, savenum)

cell_types = {'stellate','pyramidal','fast spiking'};
% experiments = {'excitation','inhibition'};
experiments = {'inhibition'};

% dataTypes = {'power','frequency','phase'};
dataTypes = {'totalGammaPower','totalFastPower'};
% dataTypes = {'frequencyACF'};

% for total power
gammaRange = [60 140];
fastRange = [100 200];

wavelet_norm = 'L2'; % L1 (default) or L2

% find ID parameters for analysis
params.locations = 'all';
params.cell_nums = 'all';
params.protocols = {'theta'};
protocol = params.protocols{1};

if strcmp(dataSets,'all'); dataSets = {'Camk2','Thy1','PV Transgenic','PV Viral'}; end

nCommentGroups = 2;
nCellTypes = length(cell_types);
nExperiments = length(experiments);
nDataTypes = length(dataTypes);

for iSet = 1:numel(dataSets)
dataSet = dataSets{iSet};

[info,~,data_path] = getInfo(dataSet);

switch dataSet
    case 'Thy1'
         saveFolder = 'C:\Users\brndn\Downloads\Thy1-ChR2\Raw Data\mEC\results\Summary';
%          saveFolder = 'C:\Users\brndn\Downloads\Thy1-ChR2\Raw Data\mEC\results\Summary\FS Only';
    case 'PV Transgenic'
         saveFolder = 'C:\Users\brndn\Downloads\PV-ChR2 Transgenic\Summary';
    case 'PV Viral'
         saveFolder = 'C:\Users\brndn\Downloads\PV-ChR2\Summary';
    case 'Camk2'
         saveFolder = 'C:\Users\brndn\Downloads\CaMK2-ChR2\Summary';
end

if ~exist(saveFolder, 'dir')
   mkdir(saveFolder)
end


for iExp = 1:nExperiments
params.experiments = experiments(iExp);
experiment = experiments{iExp};

for iType = 1:nDataTypes
dataType = dataTypes{iType};

allData = cell(nCellTypes,nCommentGroups);
diffData = cell(nCellTypes,1);
percentDiff = cell(nCellTypes,1);

for iCell = 1:nCellTypes
    cellType = cell_types{iCell};
    params.cell_types = cell_types(iCell);
    comments = cell(nCommentGroups,1);

    for iComment = 1:nCommentGroups

        if iComment == 1
            params.comments = comment_group1;
        elseif iComment == 2
            params.comments = comment_group2;
        end
        comments{iComment} = cell2mat(params.comments);

        IDs = getIDs(info,params);

        IDs = removeIDs(IDs,info,dataSet); % skip bad recordings - get params for cells to skip
        
        if (isempty(IDs)); fprintf('No %s %s Files Found.',params.cell_types,comments{iComment}); continue; end
        
        nIDs = numel(IDs);

        dataArray = zeros(nIDs,1);

        for iID = 1:nIDs
        
            ID = IDs{iID};
            ID_index = find_index(info,'ID',ID);

            p = info(ID_index);
        
            if isempty(p.comments)
                comment = 'No Comment';
            else
                comment = p.comments;
            end

            filename = sprintf('%s.abf',ID);
            fprintf('%s %s %s Analyzing %s,File %d/%d\n',dataSet,cellType,comment,filename,iID,nIDs)
        
            commentFolder = sprintf('%sresults\\%s\\%s\\%s\\%s\\%s\\%s\\%s',data_path,p.location,p.cell_type,p.cell_num,p.experiment,p.protocol,comment);
        
            dataFolder = getIDFolder(commentFolder,ID);
        
            dataFilename = [dataFolder filesep 'data.mat'];
        
            load(dataFilename,'file');
    
            % gather relevant data 
%             phase = file.cell.CWTmaxValues(1);
%             freq = file.cell.CWTmaxValues(2);
%             power = file.cell.CWTmaxValues(3);

            if any(strcmpi(dataType,{'totalGammaPower','totalFastPower'}))
                % total power
                fieldName = 'cell';
                [x,y,z,file.(fieldName).CWTcycleValues,coi] = thetaCWT(file,fieldName, file.cycle_start_index_noArtifacts,wavelet_norm);
                file.(fieldName).CWTmaxValues = getMaxValues(x,y,z);
                file.(fieldName).CWTpeakStats = getXmaxPeakStats(x,y,z,file.(fieldName).CWTmaxValues(1),[0 200],0);
                file.(fieldName).CWTsnr = file.(fieldName).CWTmaxValues(3)/mean(z(:)); % peak power / average power
        
                % only sum over regions in the cone of influence (COI) - ignores boundary effects of wavelet transform
                coi = reshape(coi,1,[]); % (1 x N) % coi represents the lowest valid frequency for each time point
                coiMask = y >= coi;   % results in [numFreq x N] logical mask
                z = z.*coiMask; % zero wavelet coefficients outside COI
        
                f = y(:,1); % frequency axis of scalogram (nFreqs x 1)
            end

%             % sum power over different frequency ranges
%             gammaMask = (f >= gammaRange(1) & f <= gammaRange(2));
%             totalGammaPower = sum(sum(z(gammaMask,:)));  % integrate over gamma frequencies and time
%     
%             fastMask = (f >= fastRange(1) & f <= fastRange(2));
%             totalFastPower = sum(sum(z(fastMask,:)));  % integrate over gamma frequencies and time

            switch dataType
                case 'power'
                    power = file.cell.CWTmaxValues(3);
                    data = power; 
                case 'frequency'
                    freq = file.cell.CWTmaxValues(2);
                    data = freq;
                case 'phase'
                    phase = file.cell.CWTmaxValues(1);
                    data = phase;
                case 'totalGammaPower'
                    % sum power over different frequency ranges
                    gammaMask = (f >= gammaRange(1) & f <= gammaRange(2));
                    totalGammaPower = sum(sum(z(gammaMask,:)));  % integrate over gamma frequencies and time
                    data = totalGammaPower;
                case 'totalFastPower'
                    fastMask = (f >= fastRange(1) & f <= fastRange(2));
                    totalFastPower = sum(sum(z(fastMask,:)));  % integrate over gamma frequencies and time
                    data = totalFastPower;
                case 'frequencyACF'
                    % --- Estimate with autocorrelation
                    plotnum = 0;
                    cycles = 2:11;
                    nCycles = numel(cycles);
                    acfFreqs = zeros(nCycles,1);
                    peakValues = acfFreqs;
                    welchValues = acfFreqs;
                    for iCycle = 1:nCycles
                       cycle = cycles(iCycle);
                       cycle_start = file.cycle_start_index_noArtifacts(cycle);
                       cycle_stop = cycle_start + file.cycle_length; 
                       
                       window_data = file.cell.gamma_data(cycle_start:cycle_stop);
                        [acfFreqs(iCycle), ~, infoACF,~] = acf_freq_estimate(window_data, file.Fs, plotnum, ...
                            'FreqRange', [50 200], ...        % Hz (set empty [] to disable)
                            'Bandpass',  [], ...        % Hz (set [] to skip)
                            'Detrend',   false, ...
                            'MaxLag',    0.02, ...            % seconds of ACF to compute
                            'MinPeakProminence', 0); % corr: first positive peak - first negative peak
                        peakValues(iCycle) = infoACF.peakValue;
                        welchValues(iCycle) = infoACF.fftpeak;
%                         fprintf('ACF estimate: f0 = %.2f Hz, peak r = %.3f\n', acfFreqs(iCycle), peakValues(iCycle));
                    end
                    medACF = median(acfFreqs,'omitnan');
                    meanR = mean(peakValues,'omitnan');
                    medFFT = median(welchValues,'omitnan');
                    fprintf('ACF estimate: f0 = %.2f Hz, peak r = %.3f, Welch estimate: f = %.2f Hz\n', medACF, meanR, medFFT);
                    data = medACF;
%                     data = medFFT;
            end

            if ~isempty(data)
                dataArray(iID,:) = data;
            end

            clearvars file;
        end

        allData{iCell,iComment} = dataArray; % save cell data
    end

end

for idiff = 1:nCellTypes
    percentDiff{idiff} = 100*((allData{idiff,2} - allData{idiff,1}) ./ allData{idiff,1}); % norm
end

% log transform allData for plotting and diffData for stats and plotting
if any(strcmp(dataType,{'power','totalGammaPower','totalFastPower'})) % convert to log scale
    allData = cellfun(@log10,allData,'UniformOutput',false); % output as cell array
end

% save before and after data in separate files
labels = {'DNQX before','DNQX after'};
for iComment = 1:nCommentGroups
    commentData = allData(:,iComment);
    data_table = myDataTable(commentData,cell_types);
    summary = mySummaryTable(commentData,cell_types); % calculate mean, SEM and N 
    sw = mySWstats(commentData,cell_types,comments);
    saveFilename = [saveFolder filesep sprintf('%s %s %s %s %s data.xlsx',dataSet,protocol,experiment,dataType,labels{iComment})];
    writetable(data_table,saveFilename,'Sheet','Data','WriteRowNames',true,'WriteMode','overwritesheet');  % save data table
    writetable(summary,saveFilename,'Sheet','Summary','WriteMode','overwritesheet','WriteRowNames',true);  % save summary stats table
    writetable(sw.results,saveFilename,'Sheet','Shapiro-Wilk','WriteRowNames',true,'WriteMode','overwritesheet');  % save sw stats table
end

for idiff = 1:nCellTypes
    diffData{idiff} = allData{idiff,2} - allData{idiff,1}; % diff of log data for plotting - stats are nonparametric
end

% dependent group comparisons
saveFilenamePaired = [saveFolder filesep sprintf('%s %s %s %s DNQX Paired',dataSet,protocol,experiment,dataType)];
statsPaired = myPairedStats(diffData,cell_types,{'Paired DNQX'},[saveFilenamePaired '.xlsx'],savenum);
% statsPaired = myPairedStats(percentDiff,cell_types,{'Paired DNQX'},[saveFilenamePaired '.xlsx'],savenum);

% independent group comparisons - use non parametric stats for differential data
saveFilenameUnpaired = [saveFolder filesep sprintf('%s %s %s %s DNQX Unpaired',dataSet,protocol,experiment,dataType)];
statsUnpaired = myMultipleIndependentGroupStats(diffData,cell_types,comments,'linear',[saveFilenameUnpaired '.xlsx'],savenum);
% statsUnpaired = myMultipleIndependentGroupStats(percentDiff,cell_types,comments,[saveFilenameUnpaired '.xlsx'],savenum);

% % log transform allData for plotting and diffData for stats and plotting
% if strcmp(dataType,'power') % convert to log scale
%     allData = cellfun(@log10,allData,'UniformOutput',false); % output as cell array
% end

% paired plots
pairedTitle = sprintf('%s %s %s DNQX Paired',dataSet,protocol,experiment);
groupLabels = myPairedGroupLabels(cell_types);
[colors,alphas] = getColorAlpha(cell_types,nCommentGroups);
y = getYparams(dataType,protocol,'normal');

figPaired = myPairedBoxplot(allData,groupLabels,y,colors,alphas,statsPaired.t2.p);
sgtitle(pairedTitle,'Fontweight','bold');

if savenum
    print(figPaired,'-vector','-dsvg',[saveFilenamePaired '.svg']); % save plot as svg file
end
 
% % log transform diffData for plotting
% if strcmp(dataType,'power') % convert to log scale
%     diffData = cellfun(@(x)sign(x).*log10(abs(x)),diffData,'UniformOutput',false); % output as cell array
% end

% independent comparisons
unpairedTitle = sprintf('%s %s %s DNQX Unpaired',dataSet,protocol,experiment);
groupLabels = myGroupLabels(cell_types);
[colors,alphas] = getColorAlpha(cell_types,1);
ydiff = getYparams(dataType,protocol,'difference');

statsUnpaired.kw.sig = 'symbol';
figBox = myBoxplot(diffData,groupLabels,ydiff,colors,alphas,statsUnpaired.kw);
sgtitle(unpairedTitle,'Fontweight','bold');

if strcmp(dataType,'power') || strcmp(dataType,'frequency')
    figBar = myBarChart(diffData,groupLabels,ydiff,colors,alphas,statsUnpaired.kw);
    sgtitle(unpairedTitle,'Fontweight','bold');
    yNorm = getYparams(dataType,protocol,'percentdiff');
    figBarNorm = myBarChart(percentDiff,groupLabels,yNorm,colors,alphas,statsUnpaired.kw);
    sgtitle([unpairedTitle ' percentdiff'],'Fontweight','bold');
end

figViolin = plotViolin(diffData,groupLabels,ydiff,colors,alphas,statsUnpaired.kw);
sgtitle(unpairedTitle,'Fontweight','bold');

if savenum
    print(figBox,'-vector','-dsvg',[saveFilenameUnpaired ' Boxplot.svg']); % save boxplot as svg file
    if strcmp(dataType,'power') || strcmp(dataType,'frequency')
        print(figBar,'-vector','-dsvg',[saveFilenameUnpaired ' BarChart.svg']); % save barchart as svg file
        print(figBarNorm,'-vector','-dsvg',[saveFilenameUnpaired ' PercentDiff BarChart.svg']); % save barchart as svg file
    end
    print(figViolin,'-vector','-dsvg',[saveFilenameUnpaired ' Violin.svg']); % save violin plots as svg file
end



end

end

end



function stats = myPairedStats(allData,rowLabels,colLabels,saveFilename,savenum)

stats.data_table = myDataTable(allData,rowLabels);
stats.summary = mySummaryTable(allData,rowLabels); % calculate mean, SEM and N 
stats.sw = mySWstats(allData,rowLabels,colLabels);
stats.levene = myVarTest(allData,rowLabels,colLabels);

stats.t2 = myPairedTtestStats(allData,rowLabels);
stats.sr = myPairedSignRankStats(allData,rowLabels);

if savenum
    if ~isempty(stats)
        writetable(stats.data_table,saveFilename,'Sheet','Data','WriteRowNames',true,'WriteMode','overwritesheet');  % save data table
        writetable(stats.summary,saveFilename,'Sheet','Summary','WriteMode','overwritesheet','WriteRowNames',true);  % save summary stats table
        writetable(stats.t2.results,saveFilename,'Sheet','Paired T-test','WriteRowNames',true,'WriteMode','overwritesheet');  % save t2 stats table
        writetable(stats.sr.results,saveFilename,'Sheet','Sign Rank test','WriteRowNames',true,'WriteMode','overwritesheet');  % save sr stats table
        writetable(stats.sw.results,saveFilename,'Sheet','Shapiro-Wilk','WriteRowNames',true,'WriteMode','overwritesheet');  % save sw stats table
        writetable(stats.levene.results,saveFilename,'Sheet','Levene','WriteRowNames',true,'WriteMode','overwritesheet');  % save levene stats table
    end
end

end



function t2 = myPairedTtestStats(data_all,cell_types)

cell_types = reshape(cell_types,[],1); % column vector

nCells = numel(cell_types);
t2.p = zeros(nCells,1); % ttest
t2.tstat = zeros(nCells,1); % ttest
t2.df = zeros(nCells,1); % ttest

for i = 1:nCells
    [~,t2.p(i),~,stats] = ttest(data_all{i},0,'tail','both');
    t2.tstat(i) = stats.tstat;
    t2.df(i) = stats.df;
end
t2.results = array2table([t2.p t2.tstat t2.df]); % set up table 
t2.results = addvars(t2.results,cell_types,'Before',1);
t2.results.Properties.VariableNames = {'Cell Types','P-value','T-statistic','df'};

end

function sr = myPairedSignRankStats(data_all,cell_types)

cell_types = reshape(cell_types,[],1); % column vector

nCells = numel(cell_types);
sr.p = zeros(nCells,1); % signrank
sr.W = zeros(nCells,1); % signrank

for i = 1:nCells
    [sr.p(i),~,stats] = signrank(data_all{i},0,'tail','both','method','exact');
    sr.W(i) = stats.signedrank; 
end
sr.results = array2table([sr.p sr.W]); % set up table
sr.results = addvars(sr.results,cell_types,'Before',1);
sr.results.Properties.VariableNames = {'Cell Types','P-value','W'};

end


function groupLabels = myPairedGroupLabels(cell_types)

nCells = numel(cell_types);

groupLabels = cell(2*nCells,1);
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
    groupLabels(2*i-1) = {sprintf('%s_{CTRL}',name)};
    groupLabels(2*i)  = {sprintf('%s_{DNQX}',name)};
%     groupLabels(2*i-1) = {'CTRL'};
%     groupLabels(2*i)  = {'DNQX'};
end


end


function fig = myPairedBoxplot(data_all,groupLabels,y,colors,alphas,t2_p)

colors = flip(colors);
alphas = flip(alphas);

tickfontsize = 15;
pvaluefontsize = 20;

[nCells,nComms] = size(data_all);

switch nCells
    case 1
        xPos = [0.95 1.05];
    case 2
        xPos = [0.95 1.05 1.25 1.35];
    case 3
        xPos = [0.95 1.05 1.25 1.35 1.55 1.65];
end

data_trans = data_all'; % for easier grouping transpose so groups are in order
maxNumEl = max(cellfun(@numel,data_trans(:)));
data_all_pad = cellfun(@(x){padarray(x(:),[maxNumEl-numel(x),0],NaN,'post')}, data_trans(:)); % Pad each vector with NaN values to equate lengths
data_all_mat = cell2mat(data_all_pad); 

data_groups = zeros(numel(data_all_mat),1);
count = 0;
for ii = 1:nCells
    for iii = 1:nComms
        count = count + 1;
        ix1 = (count-1).*(maxNumEl) + 1;
        ix2 = ix1 + maxNumEl - 1;
        data_groups(ix1:ix2) = count.*ones(maxNumEl,1);
    end
end


fig = figure;

if ~isempty(data_all_mat)
    plot_handle = boxplot(data_all_mat,data_groups,'Colors','k','Symbol','k.','OutlierSize',10,'Positions',xPos);
    set(findobj(gcf,'-regexp','Tag','\w*Whisker'),'LineStyle','-')
    set(plot_handle, 'Linewidth',1.5)
    xticks(xPos)
    xticklabels(groupLabels)
    h = findobj(gca,'Tag','Box');
    for j=1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),colors{j},'FaceAlpha',alphas{j});
    end
    drawnow;
    set(findobj(gcf, 'type', 'line', 'Tag', 'Median'), 'Color', 'k');    
end
box off

hold on

% plot data markers
% count = 0;
% for ii = 1:nComments
%     for iii = 1:nCellTypes
%           count = count + 1;
%           plot(boxplot_positions(count)*ones(length(data_all_pad{count}),1),data_all_pad{count},'xk','Markersize',5)
%     end
% end

% plot lines connecting data pairs

for iii = 1:nCells
    line_index = nComms*(iii-1) + (1:nComms);
    line_comment = xPos(line_index);
    data_before = data_all{iii,1};
    data_after = data_all{iii,2};
    for iData = 1:numel(data_before)
        line_data = [data_before(iData),data_after(iData)];
        plot(line_comment,line_data,'-','Color',0.3*ones(3,1),'Linewidth',0.8)
    end
end
hold off

ax = gca(); 
ax.TickLabelInterpreter = 'tex';  % needed for some plots like boxplot.

ax.XAxis.FontSize = tickfontsize;
ax.XAxis.FontWeight = 'bold';
xtickangle(ax,45)

ylabel(y.labelstring,'FontSize',25,'FontWeight','bold')
ax.YLim = [y.min y.max];
ax.YAxis.Scale = y.scale; 
ax.YTick = y.ticks;
ax.YTickLabel = y.tickLabels;
ax.YAxis.FontSize = tickfontsize;
ax.YAxis.FontWeight = 'bold';

% plot stats

hold on
for i = 1:nCells
    if t2_p(i) < 0.05
        data_max = max(data_all_mat(:));
        if strcmp(ax.YAxis.Scale,'log')
            horz_line_height = 10^(0.90*( log10(y.max)-log10(y.min) ) + log10(y.min));
            vert_line_min = 10^(0.87*( log10(y.max)-log10(y.min) ) + log10(y.min));
            if  data_max > vert_line_min
                vert_line_min = 10^(0.1*( log10(y.max)-log10(data_max) ) + log10(data_max));
                horz_line_height = 10^(0.2*( log10(y.max)-log10(data_max) ) + log10(data_max));
            end
            log_array = logspace(log10(horz_line_height),log10(y.max),3);
            text_y = log_array(2); % get log spaced middle height
        else
            horz_line_height = y.min + 0.90*(y.max-y.min);
            vert_line_min = y.min + 0.87*(y.max-y.min);
            if  data_max > vert_line_min
                vert_line_min = data_max + 0.1*(y.max - data_max);
                horz_line_height = data_max + 0.2*(y.max - data_max);
            end
            text_y = horz_line_height + 0.3*(y.max - horz_line_height);
        end
    
        x1 = xPos(2*i-1);
        x2 = xPos(2*i);
    
        plot( [x1 x2],horz_line_height*ones(1,2),'-k','Linewidth',1.2 )
        plot( [x1 x1],[vert_line_min horz_line_height],'-k','Linewidth',1.2 )
        plot( [x2 x2],[vert_line_min horz_line_height],'-k','Linewidth',1.2 )
        
    
        if t2_p(i) < 0.05 && t2_p(i) > 0.01
            text_x = x1 + 0.4*(x2-x1);
            text_box = text( text_x, text_y, '*','Fontweight','bold');
        elseif t2_p(i) < 0.01 && t2_p(i) > 0.001
            text_x = x1 + 0.3*(x2-x1);
            text_box = text( text_x, text_y, '**','Fontweight','bold' );
        else 
            text_x = x1 + 0.2*(x2-x1);
            text_box = text( text_x, text_y, '***','Fontweight','bold');
        end
        text_box.FontSize = pvaluefontsize;
    end
end
hold off


end

function data_table = myDataTable2(data_all,cell_types)

    nCells = numel(cell_types);
    
    % construct data array with NaNs
    maxNumEl = max(cellfun(@numel,data_all(:)));
    data_all_pad = cellfun(@(x){padarray(x(:),[maxNumEl-numel(x),0],NaN,'post')}, data_all(:)); % Pad each vector with NaN values to equate lengths
    data_all_mat = cell2mat(data_all_pad'); 
    
    % build data table
    data_table = table;

    for i = 1:nCells
        switch cell_types{i}
            case {'stellate','Stellate'}
                name = 'Stellate';
            case {'pyramidal','Pyramidal'}
                name = 'Pyramidal';
            case {'fast spiking','FastSpiking'}
                name = 'FastSpiking';
            otherwise
                name = cell_types{i};
        end
        data_table.(name) = data_all_mat(:,i);
    end

    end


end

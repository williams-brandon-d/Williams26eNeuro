function plotPatchThetaSpikesPolar(dataSets,dataTypes, savenum)

titlefontsize = 15;
cell_types = {'stellate','pyramidal','fast spiking'};

% find ID parameters for analysis
params.locations = 'all';
params.experiments = {'currentclamp'};
params.cell_nums = 'all';
params.comments = {'','DNQX before','GBZ before' 'DNQX_GABAzine_before', 'Gabazine before'}; % 'GBZ before' 'DNQX_GABAzine_before', 'Gabazine before', 'DNQX before'
params.protocols = {'theta','theta_2chan'};

if strcmp(dataSets,'all'); dataSets = {'Camk2','Thy1','PV Transgenic','PV Viral'}; end

% load data for figure

nCellTypes = numel(cell_types);
nDataTypes = numel(dataTypes);

for iSet = 1:numel(dataSets)
dataSet = dataSets{iSet};

[info,~,data_path] = getInfo(dataSet);

switch dataSet
    case 'Thy1'
         saveFolder = 'C:\Users\brndn\Downloads\Thy1-ChR2\Raw Data\mEC\results\Summary';
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

gammaCell = cell(nCellTypes,1);

for iCell = 1:nCellTypes
    cellType = cell_types{iCell};
    params.cell_types = cell_types(iCell);
    IDs = getIDs(info,params);
    
    if (isempty(IDs)); fprintf('No %s Files Found.',params.cell_types); continue; end

    validIDs = removeBadLFP(IDs,info,dataSet); % remove bad LFP recordings (CamK2 only)
    
    nIDs = numel(IDs);

    d = struct;
    
    for iID = 1:nIDs
    
        ID = IDs{iID};
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
        phases = file.theta_spike_phases; % phases
        FRs = file.FRs;
        cycle_idx_cell = file.theta_spike_indices;

        removeIndices = false(numel(FRs),1);

        removeNum = 2; % num of beginning cycles to remove from analysis
        if removeNum
            removeIndices(1:removeNum) = 1; % add to removeIndices
        end
        
        % remove empty cells
        emptyIdx = cellfun('isempty',phases); % empty idx
        removeIndices = removeIndices | emptyIdx; % add empty cells to removeIndices

        % remove burst spikes for PPC analysis - keep first spike in burst
        FR_threshold = 150; % (Hz) firing rate threshold for burst firing
        for i = 1:numel(phases)
            keepIdx = find(FRs{i} < FR_threshold) + 1;
            keepIdx = reshape(keepIdx,1,[]);
            cycle_phases = phases{i};
            if isempty(cycle_phases); continue; end
            phases{i} = cycle_phases([1 keepIdx]); % keep consecutive spikes < threshold rate
            cycle_idx = cycle_idx_cell{i};
            cycle_idx_cell{i} = cycle_idx([1 keepIdx]);  
        end

        for iType = 1:nDataTypes % stim, theta, gamma
            dataType = dataTypes{iType};

            switch dataType
                case 'stim'
                    phases = cellfun(@(x) reshape(x,[],1),phases(~removeIndices),'UniformOutput',false); % reshape all arrays in cells to columns
                    d.(dataType).phases{iID} = vertcat(phases{:}); % concatenate column data
                    
                case {'theta','gamma'}
                    if strcmp(p.protocol,'theta_2chan') % only 2 channel contains LFP
                        % check if bad lfp data
                        if any(strcmp(ID,validIDs))
                            phases = getSpikePhase(file,cycle_idx_cell,removeIndices,dataType);
                            phases = cellfun(@(x) reshape(x,[],1),phases,'UniformOutput',false); % reshape all arrays in cells to columns
                            d.(dataType).phases{iID} = vertcat(phases{:}); % concatenate column data
                        else
                            d.(dataType).phases{iID} = {};
                        end
                    else
                        d.(dataType).phases{iID} = {};
                    end
            end

        end
    
        clearvars file;
    end

    for iType = 1:nDataTypes
        dataType = dataTypes{iType};
        dataCell = d.(dataType).phases;

        dataCell(cellfun(@isempty,dataCell)) = []; % remove empty cells

        % save gamma phases for each cell type
        if strcmp(dataType,'gamma')
            gammaCell{iCell} = dataCell;
        end

        nNeurons = numel(dataCell);

        [color,~] = getColorAlpha({cellType},1);
        color = color{1};

        figCell = plotPolar(dataCell,color); % plot figure for each cell type
        ax = gca;
        switch dataType
            case {'stim','theta'}
                ax.RLim = [0 0.20];
                ax.RTick = [0.1 0.2];
            case 'gamma'
                ax.RLim = [0 0.10];
                ax.RTick = [0.05 0.1];
        end
        titleString = sprintf('%s %s %s Spike Phase (n=%d)',dataSet,cellType, dataType, nNeurons);
        title(titleString,'FontSize',titlefontsize,'Fontweight','bold');
        
        if savenum
            saveFilename = [saveFolder filesep sprintf('%s %s %s Spike Phase.svg',dataSet,cellType,dataType)];
            print(figCell,'-vector','-dsvg',saveFilename);
        end
        
    end

end

% plot resultant vector for gamma all cell types
VS = NaN(nCellTypes,1);
mu = VS;
Z = VS;
colors = NaN(nCellTypes,3);

for iCell = 1:nCellTypes
    phaseCell = gammaCell{iCell};
    phases_cell = cellfun(@(x) reshape(x,[],1),phaseCell,'UniformOutput',false); % reshape all arrays in cells to columns
    phases_mat = vertcat(phases_cell{:}); % concatenate column data
    
    % calculate vector strength and angle
    [VS(iCell), mu(iCell), Z(iCell)] = vector_strength(phases_mat);

    [color,~] = getColorAlpha(cell_types(iCell),1);
    colors(iCell,:) = color{1};    
end

% 0 degrees at the top and rotate clockwise like polar plots
Z  = 1i * conj(Z);

% plot 
figVS = figure;
c = compass(Z);
hold on;
for iCell = 1:nCellTypes
    c1 = c(iCell);
    c1.LineWidth = 2;
    c1.Color = colors(iCell,:);
end

compass_add_cardinal_rad_labels(gca, 'pi');  % 'pi' or 'dec'

hold off;
% titleString = sprintf('%s Gamma Vector Strength',dataSet);
% title(titleString,'FontSize',titlefontsize,'Fontweight','bold');

% build table of VS and angles

% column names for cell types
varNames = cell(nCellTypes,1);
for i = 1:nCellTypes
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
    varNames{i} = name;
end

data_table = array2table([VS'; mu';],'VariableNames',varNames,'RowNames',{'VS','Phase'});

if savenum
    saveFilename = [saveFolder filesep sprintf('%s Gamma Vector Strength',dataSet)];
    print(figVS,'-vector','-dsvg',saveFilename);
    writetable(data_table,[saveFilename '.xlsx'],'WriteMode','overwrite','WriteRowNames',true);  % save summary stats table
end


end


function phases = getSpikePhase(file,cycle_idx_cell,removeIndices,dataType)
    % load spike indices and lfp phase
%     cycle_idx_cell = file.theta_spike_indices;
    cycle_idx_cell = cycle_idx_cell(~removeIndices);

    switch dataType
        case 'theta'
            lfp_phase = file.lfp.theta_phase;
        case 'gamma'
            lfp_phase = file.lfp.gamma_phase;
    end

    % get spike-lfp phases
    cycle_start_indices = file.cycle_start_index_noArtifacts(~removeIndices(file.cycles_noArtifacts));

    % get spike-lfp phases
    nCycles = numel(cycle_start_indices);
    phases = cell(nCycles,1);

    for icycle = 1:nCycles
        cycle_indices = cycle_start_indices(icycle) + (0:file.cycle_length);
        cycle_phase = lfp_phase(cycle_indices);
        cycle_spike_indices = cycle_idx_cell{icycle,:};
        phases{icycle,:} = cycle_phase(cycle_spike_indices)';
    end
    

end


end

function compass_add_cardinal_rad_labels(ax, mode)
% Add cardinal (every 90°) tick labels in radians for COMPASS with
% 0 at top and clockwise increase. Works even if no degree labels exist.
%   compass_add_cardinal_rad_labels(gca)        % π-style labels
%   compass_add_cardinal_rad_labels(gca,'dec')  % decimal radians

    if nargin < 1 || isempty(ax), ax = gca; end
    if nargin < 2, mode = 'pi'; end

    % Where to place labels: just inside/outside the circle
    % (auto-scales if you didn't use a unit dummy)
    r = max([abs(xlim(ax)) abs(ylim(ax))]);
    rlab = 1.07 * r;   % push slightly outside the circle

    % Cardinal angles in north-CW radians
    th = [0, pi/2, pi, 3*pi/2];                     % 0=N, π/2=E, π=S, 3π/2=W
    alpha = pi/2 - th;                               % display angle used by COMPASS
    x = rlab * cos(alpha);  y = rlab * sin(alpha);

    % Labels (π style or decimals)
    if strcmpi(mode,'dec')
        lbls = {'0','1.57','3.14','4.71'};
        interp = 'tex';
    else
        lbls = {'0','π/2','π/-π','-π/2'};        % TeX π
        interp = 'none';
    end

    % remove integer labels
    txt = findall(ax,'Type','text');
    for k = 1:numel(txt)
        s = string(txt(k).String);
        if ~isempty(regexp(s, '^\s*\d+\s*(?:°|deg)?\s*$', 'once'))
            delete(txt(k));
        else % make bold

        end
    end

    % Draw new labels
    for k = 1:numel(th)
        text(x(k), y(k), lbls{k}, 'Parent', ax, ...
            'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
            'Interpreter', interp, 'Clipping','off', 'Tag','compassCardinal','FontSize',25,'FontWeight','bold');
    end
end
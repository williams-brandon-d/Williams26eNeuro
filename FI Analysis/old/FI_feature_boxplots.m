%% compare features across excitatory cells
close all
clear variables
clc

% stellate id 21819031 has long latency
% stellate id 20901000 and 20910000 have long time constant

cd('C:\Users\brndn\OneDrive\Desktop\White Lab\thy1_chr2\Thy1-ChR2\MATLAB\thy1chr2\')

info_path = 'C:\Users\brndn\OneDrive\Desktop\White Lab\thy1_chr2\Thy1-ChR2\MATLAB\thy1chr2\';
% info_name = 'thy1_chr2_info_082021.mat';
info_name = 'thy1_chr2_info_032822.mat';

load(fullfile(info_path,info_name),'info')

p.experiment = 'currentclamp';
p.protocol = 'FI';

p.location = 'all';
cell_types = {'stellate','pyramidal'};
p.cell_num = 'all';
comments = {'','DNQX before'}; % combine comment data

% features = {'dap'};
% plotTitles = {'depolarizing After Potential (mV)'};
features = {'ISI_ratio','sag_hyper','sag_depol_thresh','dap','latency','time_constant'};
plotTitles = {'ISI Ratio 1/2';
              'Sag hyperpolarizing (mV)';
              'Sag depolarizing (mV)'; 
              'depolarizing After Potential (mV)';
              'Latency to first spike (ms)';
              'Membrane Time Constant (ms)'};
tickfontsize = 15;
titlefontsize = 20;

nFeatures = length(features);
nTypes = length(cell_types);
nComments = length(comments);

feature = cell(nTypes,1);
row1 = feature';
row2 = row1;
tickLabels = feature;

ks_h = zeros(nFeatures,nTypes);
ks_p = ks_h;
means = ks_h;
stds = ks_h;

t2_h = zeros(nFeatures,1);
t2_p = t2_h;

IDs_comments = cell(nComments,1);

for iFeature = 1:nFeatures
    for iType = 1:nTypes
        p.cell_type = cell_types{iType};
        
        for j = 1:length(comments)
            p.comments = comments{j};
            IDs_comments{j} = find_IDs(info,p);
        end
        IDs = vertcat(IDs_comments{:});
        nIDs = length(IDs);
        ID_index = zeros(nIDs,1);
        feature_array = zeros(nIDs,1);
        
        for m = 1:nIDs
            ID_index(m) = find_index(info,'ID',IDs{m});
            info_ID = info(ID_index(m));
            get_feature = eval( sprintf( 'info_ID.%s',features{iFeature}) );
            if isempty(get_feature)
                feature_array(m) = NaN;
            else
                feature_array(m) = get_feature;
            end
        end
        feature{iType} = feature_array;
        
        means(iFeature,iType) = mean(feature{iType});
        stds(iFeature,iType) = std(feature{iType});

        if numel(feature{iType}) > 1 && ((stds(iFeature,iType)) > 0)
            normData = ( feature{iType} - means(iFeature,iType) )/stds(iFeature,iType);
            [ks_h(iFeature,iType),ks_p(iFeature,iType)] = kstest(normData);  % 0 = normal, 1 = not
        end
        
        row1{iType} = sprintf('%s',p.cell_type);
        row2{iType} = sprintf('(n=%d)',length(feature{iType}));
    end
    
    labelArray = [row1; row2]; 
    labelArray = strjust(pad(labelArray),'center'); % 'left'(default)|'right'|'center'
    tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
    
    [t2_h(iFeature),t2_p(iFeature)] = ttest2(feature{1},feature{2}); % compare feature between cell types
    
    maxNumEl = max(cellfun(@numel,feature));
    feature_pad = cellfun(@(x){padarray(x(:),[maxNumEl-numel(x),0],NaN,'post')}, feature); % Pad each vector with NaN values to equate lengths
    feature_mat = cell2mat(feature_pad); 
    feature_mat = reshape(feature_mat,[],nTypes);
    
    figure;
    h = boxplot(feature_mat);
    set(h, 'Linewidth',1.5)
    hold on
    for plotType = 1:nTypes
        plot(plotType*ones(length(feature{plotType}),1),feature{plotType},'xk')
    end
    ylimits = ylim;
    ymax = ylimits(2);
    plot( [1 2],1.1*[ymax ymax],'-k','Linewidth',1.2 )
    plot( [1 1],[1.05*ymax 1.1*ymax],'-k','Linewidth',1.2 )
    plot( [2 2],[1.05*ymax 1.1*ymax],'-k','Linewidth',1.2 )
    if t2_p(iFeature) < 0.001
        p_rounded = 10^(ceil(log10( t2_p(iFeature) )));
        text( 1.33, 1.15*ymax, sprintf('p < %.1g',p_rounded) )
    else
        text( 1.33, 1.15*ymax, sprintf('p = %.3f',t2_p(iFeature)) )
    end
    title( plotTitles{iFeature},'FontSize',titlefontsize )
    ax = gca(); 
    ax.XTick = 1:nTypes; 
    ax.XLim = [0.5,nTypes+0.5];
    ax.XTickLabel = tickLabels; 
    ax.YLim = 1.2*ylim;
    ax.TickLabelInterpreter = 'tex';  % needed for some plots like boxplot.
    ax.XAxis.FontSize = tickfontsize;
    ax.XAxis.FontWeight = 'bold';
    ax.YAxis.FontSize = tickfontsize;
    ax.YAxis.FontWeight = 'bold';
   
end
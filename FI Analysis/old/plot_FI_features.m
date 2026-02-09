%% compare features across excitatory cells
close all
clear variables
clc

cd('C:\Users\brndn\OneDrive\Desktop\White Lab\thy1_chr2\Thy1-ChR2\MATLAB\thy1chr2\')

save_num = 0;
save_path = 'C:\Users\brndn\OneDrive\Desktop\White Lab\thy1_chr2\Thy1-ChR2\MATLAB\thy1chr2\Figures\03-07-21';

info_path = 'C:\Users\brndn\OneDrive\Desktop\White Lab\thy1_chr2\Thy1-ChR2\MATLAB\thy1chr2\';
info_name = 'thy1_chr2_info_032921.mat';

load(fullfile(info_path,info_name),'info')

p.location = 'all';
cell_types = {'stellate','pyramidal'};
p.experiment = 'currentclamp';
p.protocol = 'FI';
p.cell_num = 'all';
p.comments = '';
features = {'ISI_ratio','sag_hyper','sag_depol_thresh','dap','latency'};
plotTitles = {'ISI Ratio 1/2';
              'Sag hyperpolarizing (mV)';
              'Sag depolarizing (mV)'; 
              'depolarizing After Potential (mV)';
              'Latency to first spike (ms)'};
tickfontsize = 15;
titlefontsize = 20;

nFeatures = length(features);
nTypes = length(cell_types);
feature = cell(nTypes,1);
tickLabels = cell(nTypes,1);
ks_h = zeros(nFeatures,nTypes);
ks_p = ks_h;
means = ks_h;
stds = ks_h;
t2_h = zeros(nFeatures,1);
t2_p = t2_h;

for iFeature = 1:nFeatures
    for iType = 1:nTypes
        p.cell_type = char(cell_types(iType));

        IDs = find_IDs(info,p);
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
        tickLabels{iType} = p.cell_type;
        
        means(iFeature,iType) = mean(feature{iType});
        stds(iFeature,iType) = std(feature{iType});

        if numel(feature{iType}) > 1 && ((stds(iFeature,iType)) > 0)
            normData = ( feature{iType} - means(iFeature,iType) )/stds(iFeature,iType);
            [ks_h(iFeature,iType),ks_p(iFeature,iType)] = kstest(normData);  % 0 = normal, 1 = not
        end
    end
    
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
        text( 1.33, 1.15*ymax, 'p < 0.001' )
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
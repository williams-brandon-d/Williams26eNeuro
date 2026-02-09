%% compare features across excitatory cells
close all
clear variables
clc

% plot 3D features vs labels

cd('C:\Users\brndn\OneDrive\Desktop\White Lab\thy1_chr2\Thy1-ChR2\MATLAB\thy1chr2\')

info_path = 'C:\Users\brndn\OneDrive\Desktop\White Lab\thy1_chr2\Thy1-ChR2\MATLAB\thy1chr2\';

info_name = 'thy1_chr2_info_082021.mat';

load(fullfile(info_path,info_name),'info')

p.location = 'all';
cell_types = {'stellate','pyramidal'};
p.experiment = 'currentclamp';
p.protocol = 'FI';
p.cell_num = 'all';
p.comments = '';
comments = {'','DNQX before'};
features = {'ISI_ratio','sag_hyper','sag_depol_thresh','dap','latency','time_constant'};
plotTitles = {'ISI Ratio 1/2';
              'Sag hyperpolarizing (mV)';
              'Sag depolarizing (mV)'; 
              'depolarizing After Potential (mV)';
              'Latency to first spike (ms)';
              'Membrane Time Constant (ms)'};
% tickfontsize = 15;
% titlefontsize = 20;

nFeatures = length(features);
nTypes = length(cell_types);
nComments = length(comments);

feature = cell(nTypes,1);
all_features = cell(nFeatures,1);
% tickLabels = cell(nTypes,1);
IDs_comments = cell(nComments,1);
feature_cell_matrix = cell(nFeatures,nTypes);

% FI_index = find_index(info,'protocol','FI');
% stellate_index = find_index(info,'cell_type','stellate');
% pyramidal_index = find_index(info,'cell_type','pyramidal');
% FI_stellate_index = intersect(FI_index,stellate_index);
% FI_pyramidal_index = intersect(FI_index,pyramidal_index);

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
        feature_cell_matrix{iFeature,iType} = feature_array;
    end
   all_features{iFeature} = cell2mat(feature);
end

all_feature_vector = cell2mat(all_features);
nDataPoints = length(all_features{1});
feature_matrix = reshape(all_feature_vector,nDataPoints,nFeatures);

plotmatrix( feature_matrix )

ISI_ratio = feature_matrix(:,1);
sag_hyper = feature_matrix(:,2);
sag_depol = feature_matrix(:,3);
dap = feature_matrix(:,4);
latency = feature_matrix(:,5);
time_constant = feature_matrix(:,6);

color = load_fav_colors();

figure;
for iType = 1:nTypes
    h = plot3(feature_cell_matrix{1,iType},feature_cell_matrix{3,iType},feature_cell_matrix{6,iType},'.','Markersize',10);
    set(h,'Color',color(iType,:))
    set(h,'DisplayName',cell_types{iType})
    hold on;
end
hold off;
xlabel(plotTitles{1})
ylabel(plotTitles{3})
zlabel(plotTitles{6})
legend()

k = 2;
nReplicates = 100;

X = feature_matrix(:,[1 3 6]);

[ idx,centroids ] = kmeans(X,k,'Replicates',nReplicates);

plot_labels = cell(k,1);

figure;
for iCluster = 1:k
    cluster_idx = find(idx == iCluster);
    x1 = ISI_ratio( cluster_idx );
    x2 = sag_depol( cluster_idx );
    x3 = time_constant( cluster_idx );
%     Clusters(cluster_idx,:,iCluster) = [x1,x2,x3];
    plot_labels(iCluster) = {sprintf('Cluster%d',iCluster)};
    plot3(x1,x2,x3,'.','Color',color(iCluster,:),'Markersize',20)
    hold on;
end
xlabel(plotTitles{1})
ylabel(plotTitles{3})
zlabel(plotTitles{6})
legend(plot_labels)

[coeff,score,~,~,explained] = pca( X );

figure;
scatter3(score(:,1),score(:,2),score(:,3))
axis equal
xlabel('1st Principal Component')
ylabel('2nd Principal Component')
zlabel('3rd Principal Component')

figure;
scatter(score(:,1),score(:,2))
axis equal
xlabel('1st Principal Component')
ylabel('2nd Principal Component')


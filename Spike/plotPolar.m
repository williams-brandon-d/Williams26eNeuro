function fig = plotPolar(phases_cell,color)
    
    nbins = 36; % add different bins for different frequencies

    tickfontsize = 25;

    phases_cell = cellfun(@(x) reshape(x,[],1),phases_cell,'UniformOutput',false); % reshape all arrays in cells to columns
    phases_mat = vertcat(phases_cell{:}); % concatenate column data
    
    edges = linspace(-pi,pi,nbins+1);
    
    fig = figure;
    % alternatively could compute histogram previously
    polarhistogram(phases_mat,edges,'FaceColor',color,'FaceAlpha',0.6,'Normalization','probability');

    ax = gca;
    ax.ThetaZeroLocation = 'top';
    ax.ThetaDir = 'clockwise'; % 90 degrees at the right
    ax.ThetaTick = [0 90 180 270];
    ax.ThetaTickLabels = {'0','π/2','π/-π','-π/2',};
%     ax.RLim = [0 0.20];
%     ax.RTick = [0.1 0.2];
    ax.FontSize = tickfontsize;
    ax.FontWeight = 'bold';

end
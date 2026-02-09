function fig = plotScalogram(x,y,z,plotTitle,maxValues,data_units,stimType)
    fontsize = loadFontSizes();
    fontsize.tick = 20;
    fontsize.cbar = 14;
    fontweight = 'bold';

    fig = figure;

%     surf(x,y,z)

    % downsample surf plot resolution to save as vector graphics
    x_new  = linspace(-pi,pi,size(z,1));
%     x_new  = linspace(-pi,pi,100);
    z_ds = interp1(x(1,:), z', x_new, 'linear')';
    [X_new,Y_new] = meshgrid(x_new,y(:,1));

    surf(X_new,Y_new,z_ds,'EdgeColor', 'none');
    shading interp
    view(0,90)
    
    hcb = colorbar('Fontsize',fontsize.cbar,'Fontweight','bold');
%     title(hcb,sprintf('Power (%s^2)',data_units),'FontSize',fontsize.cbar,'FontWeight','bold') 

    ylim([-inf inf])
    clim([0 inf]);
%     clim([0 130.9301]) % change for dnqx after 
    colormap('hot')

    ax = gca;
    ax.XAxis.TickLabelInterpreter = 'tex';   % tex for x-axis

    switch stimType
        case 'theta'
            xlim([-inf inf])
            xticks([-pi -pi/2 0 pi/2 pi]);
            xticklabels({'-π','-π/2','0','π/2','π'});
            xLabel = 'Stim Theta Phase (rad)';
        case 'pulse'
            xlim([0 500]); % ms
            xTicks = 0:100:500;
            xticks(xTicks);
            xticklabels(string(xTicks));
            xLabel = 'Time (ms)';
        case 'noise'   
            xLabel = 'Time (ms)';
    end

    ax.YAxis.FontSize = fontsize.tick;
    ax.XAxis.FontSize = fontsize.tick;

    xlabel(xLabel,'Fontsize',fontsize.tick)
    ylabel('Frequency (Hz)','Fontsize',fontsize.tick)
    ax.YAxis.FontWeight = fontweight;
    ax.XAxis.FontWeight = fontweight;

    if nargin > 3
%         title(plotTitle,'FontSize',fontsize.title,'Interpreter','none')
    end

    if nargin > 4
        hold on
%         plot3(maxValues(1),maxValues(2),maxValues(3),'xb')
        hold off
    end

end

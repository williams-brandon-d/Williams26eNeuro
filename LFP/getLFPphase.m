function file = getLFPphase(file, plotnum, savenum, dsnum)
% get instantaneous phase of theta and gamma oscillations in LFP

% get phase of oscillation using hilbert transform
file.lfp.theta_phase = angle(hilbert(file.lfp.theta_data)); % from -pi to pi
file.lfp.gamma_phase = angle(hilbert(file.lfp.gamma_data)); % from -pi to pi

if dsnum % if downsample time is input (spikes)
    % downsample unwrapped lfp phase data using linear interpolation
    file.lfp.theta_phase_ds = wrapToPi(interp1(file.time,unwrap(file.lfp.theta_phase),file.time_ds));
    file.lfp.gamma_phase_ds = wrapToPi(interp1(file.time,unwrap(file.lfp.gamma_phase),file.time_ds));
end

% alternative method - downsample data first then hilbert transform
% lfp.theta_phase_ds = angle(hilbert(lfp.theta_data_ds)); % from -pi to pi
% lfp.gamma_phase_ds = angle(hilbert(lfp.gamma_data_ds)); % from -pi to pi

% check downsample quality
% plotds(lfp.time,lfp.theta_phase,spikes.time,lfp.theta_phase_ds);
% plotds(lfp.time,lfp.gamma_phase,spikes.time,lfp.gamma_phase_ds);

if plotnum
    thetaFig = plotHilbert(file.lfp.theta_data,file.lfp.theta_phase,file.dt);
    sgtitle('Theta','Fontweight','bold');
    gammaFig = plotHilbert(file.lfp.gamma_data,file.lfp.gamma_phase,file.dt);
    sgtitle('Gamma','Fontweight','bold');
%     plotHilbert(lfp.theta_data_ds,lfp.theta_phase_ds)
%     plotHilbert(lfp.gamma_data_ds,lfp.gamma_phase_ds)
end

if plotnum && savenum
    saveas(thetaFig,[file.saveFolder filesep 'hilbert theta.png']);
    saveas(gammaFig,[file.saveFolder filesep 'hilbert gamma.png']);
end

function fig = plotHilbert(y,phase,dt)
    % freq = instfreq(y,fs,'Method','hilbert');
    tickfontsize = 12;

    nSamples = length(y);
    dt_ms = 1000*dt; % dt in ms
    time = (0:nSamples-1)*dt_ms;

    fig = figure;

    subplot(2,1,1)
    plot(time,y)
    xlim([700 1000]);
    ylabel('uV');
    title('Filtered Data');
    ax = gca;
    ax.YAxis.FontSize = tickfontsize;
    ax.YAxis.FontWeight = 'bold';
    ax.XAxis.FontSize = tickfontsize;
    ax.XAxis.FontWeight = 'bold';

    subplot(2,1,2)
    plot(time,phase)
    xlim([700 800]);
    ylabel('Phase (rad)');
    yticks([-pi -pi/2 0 pi/2 pi]);
    yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
    title('Hilbert Phase');
    ax = gca;
    ax.YAxis.FontSize = tickfontsize;
    ax.YAxis.FontWeight = 'bold';
    ax.XAxis.FontSize = tickfontsize;
    ax.XAxis.FontWeight = 'bold';

    % subplot(3,1,2)
    % plot(freq)

    han = axes(fig,'visible','off'); 
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    han.Title.Visible='on';
    xlabel(han,'Time (ms)','FontSize',tickfontsize,'FontWeight','bold');

end

end
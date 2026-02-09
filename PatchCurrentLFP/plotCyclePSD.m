function [out,fig] = plotCyclePSD(data,Fs,cycle_start_indices,cycle_length,PSDtype,data_units,plotnum)
% plot PSDs each theta cycle of raw data

% find peak in gamma range with the greatest prominence 

gamma_range = [60 200]; % range to find gamma peaks 

% params for interpolation
flimits = [0 250]; % frequency range
df = 1; % frequency spacing
freq = flimits(1):df:flimits(2); % frequency interpolation axis

% params for mtftt
ntapers = 1; % number of tapers 

% params for welch PSD
nSegments = 1; % one window for PSD

nCycles = numel(cycle_start_indices);

% N = cycle_length + 1;
% nfft = max(2^(nextpow2(N)+params.pad),N); % all frequencies
% nfft = floor((nfft/2)) + 1; % positive frequencies

nfft = numel(freq);

S = zeros(nfft,nCycles);

for iCycle = 1:nCycles
    
    cycle_indices = cycle_start_indices(iCycle) + (0:cycle_length); % gather raw data indices
    
    data_cycle = data(cycle_indices)'; % raw data - samples x trials
    
    switch PSDtype
        case 'pwelch'
            [S(:,iCycle),freq] = welchPSD(data_cycle,Fs,nSegments,freq);
        case 'mtspec'
            [S(:,iCycle),freq] = mtspecPSD(data_cycle,Fs,ntapers,freq);
    end

end

Smean = mean(S,2);

% find gamma range
gamma_mask = (freq >= gamma_range(1) & freq <= gamma_range(2));

% Smean = 10*log10(Smean); % dB

gamma_start_index = find(gamma_mask,1);
% find most prominent gamma peak
[~, peak_indices,~,peak_proms] = findpeaks(Smean(gamma_mask));
[~,maxPromIdx] = max(peak_proms);
max_gamma_idx = peak_indices(maxPromIdx);

max_index = gamma_start_index + max_gamma_idx - 1;
max_gamma_freq = freq(max_index);
max_gamma_power = Smean(max_index);

% sum gamma power +- 15 Hz from the peak
half_bandwidth = 15; % Hz
nHalf = half_bandwidth/df; % number of samples
sum_gamma_power = sum(Smean(max_index-nHalf:max_index+nHalf));

% % find max gamma peak
% [max_gamma_peak,max_idx] = max(Smean(gamma_mask));
% max_gamma_freq = freq(gamma_start_index + max_idx - 1);
% 
% % find half max power indices left and right of max peak
% % half_max = min(Smean(gamma_mask)) + 0.5*range(Smean(gamma_mask));
% half_max = max_gamma_peak - 3; % -3 dB from peak is half power
% 
% left_mask = (freq >= gamma_range(1) & freq <= max_gamma_freq);
% left_half_index = find(Smean(left_mask) >= half_max,1);
% left_half_freq = freq( find(left_mask,1) + left_half_index - 1);
% left_half_peak = Smean( find(left_mask,1) + left_half_index - 1);
% 
% right_mask = (freq >= max_gamma_freq & freq <= gamma_range(2));
% right_half_index = find(Smean(right_mask) <= half_max,1);
% right_half_freq = freq( find(right_mask,1) + right_half_index - 2);
% right_half_peak = Smean( find(right_mask,1) + right_half_index - 2);

% fprintf('Gamma half max range: %d - %d Hz\n',left_half_freq,right_half_freq);

if plotnum
    tickfontsize = 15;

    % plot spectra
    fig = figure;
    % plot(f,S,'Color',0.5*[1 1 1])
    hold on
    plot(freq,Smean,'-k','LineWidth',1.5)
    yLim = ylim;
    % plot(left_half_freq,left_half_peak,'or',right_half_freq,right_half_peak,'ob');
    % if ~isempty(left_half_peak)
    % plot([left_half_freq left_half_freq],[yLim(1) left_half_peak],'--r');
    % end
    % if ~isempty(right_half_peak)
    % plot([right_half_freq right_half_freq],[yLim(1) right_half_peak],'--r');
    % end
    if ~isempty(max_gamma_power)
        plot([max_gamma_freq max_gamma_freq],[yLim(1) max_gamma_power],'--b');
    end

    hold off
    xlabel('Frequency (Hz)');
    ylabel(sprintf('PSD (%s^{2}/Hz)',data_units));
    % ylabel('PSD (dB)');
    xlim(gamma_range);
    % title(sprintf('Peak Gamma: %d Hz | Half Power Range: %d - %d Hz',max_gamma_freq,left_half_freq,right_half_freq),'FontSize',12);
    title(sprintf('Peak Gamma: %d Hz',max_gamma_freq),'FontSize',12);
    box off
    ax = gca;
    ax.YAxis.FontSize = tickfontsize;
    ax.YAxis.FontWeight = 'bold';
    ax.XAxis.FontSize = tickfontsize;
    ax.XAxis.FontWeight = 'bold';
else
    fig = [];
end

out.type = PSDtype;
out.S = S;
out.f = freq;
out.max_psd_gamma_freq = max_gamma_freq;
out.max_gamma_power = max_gamma_power;
out.sum_gamma_power = sum_gamma_power;

% data.max_psd_gamma_range = [left_half_freq right_half_freq];

end
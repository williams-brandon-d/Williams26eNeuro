function [data,fig] = plotallCyclePSD(data,Fs,cycle_start_indices,cycle_length)
% plot PSDs each theta cycle of raw lfp data

% find peak in gamma range with the greatest prominence 

TW = 1; % time-bandwidth product determines number of tapers

flimits = [0 200];

tickfontsize = 15;

gamma_range = [50 150];

% params for interpolation
method = 'cubic';
df = 1;
freq = flimits(1):df:flimits(2);

% params for mtftt
ntapers = 2*TW-1;
params.Fs = Fs;
params.tapers = [TW,ntapers];
params.pad = -1; % no padding
params.trialave = 0;

nCycles = numel(cycle_start_indices);

N = cycle_length + 1;
nfft = max(2^(nextpow2(N)+params.pad),N); % all frequencies
nfft = floor((nfft/2)) + 1; % positive frequencies

S = zeros(nfft,nCycles);

for iCycle = 1:nCycles
    
    % gather raw stim data indices
    stim_indices = cycle_start_indices(iCycle) + (0:cycle_length);
    
    % raw data
    data_stim = data.raw_data(stim_indices)'; % samples x trials
    
    % stim spectra
    [S(:,iCycle),f] = mtspectrumc( data_stim, params );

end

Smean = mean(S,2);

% S = 10*log10(S); % dB
Smean = 10*log10(Smean); % dB

% interpolate spectra to the new frequency axis
Smean = interp1(f, Smean, freq, method);

% find gamma range
gamma_mask = (freq >= gamma_range(1) & freq <= gamma_range(2));
gamma_start_index = find(gamma_mask,1);

% find most prominent gamma peak
[~, peak_indices,~,peak_proms] = findpeaks(Smean(gamma_mask));
[~,maxPromIdx] = max(peak_proms);
max_gamma_idx = peak_indices(maxPromIdx);
max_gamma_freq = freq(gamma_start_index + max_gamma_idx - 1);
max_gamma_power = Smean(gamma_start_index + max_gamma_idx - 1);

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
% ylabel('PSD (uV^2/Hz)');
ylabel('PSD (dB)');
xlim(flimits);
% title(sprintf('Peak Gamma: %d Hz | Half Power Range: %d - %d Hz',max_gamma_freq,left_half_freq,right_half_freq),'FontSize',12);
title(sprintf('Peak Gamma: %d Hz',max_gamma_freq),'FontSize',12);
box off
ax = gca;
ax.YAxis.FontSize = tickfontsize;
ax.YAxis.FontWeight = 'bold';
ax.XAxis.FontSize = tickfontsize;
ax.XAxis.FontWeight = 'bold';


data.max_psd_gamma_freq = max_gamma_freq;
data.S = S;
data.f = f;
% data.max_psd_gamma_range = [left_half_freq right_half_freq];

end
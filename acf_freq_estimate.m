function [f0, T0, info, fig] = acf_freq_estimate(x, fs, plotnum, varargin)
%ACF_FREQ_ESTIMATE  Frequency estimate from autocorrelation.
%   [f0, T0, info] = acf_freq_estimate(x, fs, 'Name',Value,...)
%   Inputs:
%       x  : signal (vector)
%       fs : sampling rate (Hz)
%   Name-Value (all optional):
%       'FreqRange'         : [fmin fmax] Hz to search peaks (default [])
%       'Bandpass'          : [f1 f2] Hz pre-filter via bandpass() (default [])
%       'Detrend'           : true/false (default true)
%       'MaxLag'            : max lag (seconds) for ACF (default 1.0)
%       'MinPeakProminence' : peak prominence for findpeaks (default 0.01)
%
%   Outputs:
%       f0   : estimated frequency (Hz) or NaN if none found
%       T0   : estimated period (s) corresponding to first significant peak
%       info : struct with diagnostics (acf, lags, peakValue, etc.)

    p = inputParser;
    addParameter(p, 'FreqRange', [], @(v) isempty(v) || (isnumeric(v) && numel(v)==2));
    addParameter(p, 'Bandpass',  [], @(v) isempty(v) || (isnumeric(v) && numel(v)==2));
    addParameter(p, 'Detrend',   true, @(v)islogical(v) || isnumeric(v));
    addParameter(p, 'MaxLag',    1.0, @(v)isnumeric(v) && v>0);
    addParameter(p, 'MinPeakProminence', 0.01, @(v)isnumeric(v) && v>=0);
    parse(p, varargin{:});
    opts = p.Results;

    x = x(:);
    x = x(~isnan(x));                 % drop NaNs
    if isempty(x), f0 = NaN; T0 = NaN; info = struct(); return; end

    % Optional pre-filter
    if ~isempty(opts.Bandpass)
        fbp = sort(opts.Bandpass);
        Wn = fbp / (fs/2);              % normalize cutoff frequencies
        filterOrder = 4;               % typical for oscillatory neural data
        [b, a] = butter(filterOrder, Wn, 'bandpass');
        x = filtfilt(b, a, x);         % zero-phase filtering
    end

    if opts.Detrend
        x = detrend(x, 'linear');
    end
    x = x - mean(x);

    % --- Autocorrelation (unbiased), normalized by r(0) ---
    maxLagSamp = min(numel(x)-1, round(opts.MaxLag*fs));
    [r, lags]  = xcorr(x, maxLagSamp, 'unbiased');
    % keep non-negative lags
    keep       = lags >= 0;
    r          = r(keep);
    lags       = lags(keep);
    if r(1) ~= 0
        r = r./r(1); % normalize so r(0)=1
    end

    % --- Limit search range by frequency (if provided) ---
    if ~isempty(opts.FreqRange)
        fmin = max(eps, min(opts.FreqRange));
        fmax = max(opts.FreqRange);
        lagMin = ceil(fs/fmax);   % smallest lag to consider
        lagMax = floor(fs/fmin);  % largest lag to consider
    else
        lagMin = 1; 
        lagMax = numel(lags)-1;
    end
    lagMin = max(lagMin, 1);
    lagMax = min(lagMax, numel(lags)-1);
    if lagMax <= lagMin
        f0 = NaN; T0 = NaN; info = struct('acf', r, 'lags', lags, 'peakValue', NaN); return;
    end

    rSearch   = r(lagMin:lagMax);
    lSearch   = lags(lagMin:lagMax);
    % Avoid picking very close peaks by enforcing distance ~ 1/fmax
    if ~isempty(opts.FreqRange)
        minPeakDist = round(fs / max(opts.FreqRange)); % samples
    else
        minPeakDist = 1;
    end

    [pks, locs] = findpeaks(rSearch, ...
        'MinPeakProminence', opts.MinPeakProminence, ...
        'MinPeakHeight', 0, ...
        'MinPeakDistance',   max(minPeakDist,1));

    if isempty(pks)
        f0 = NaN; T0 = NaN; 
        info = struct('acf', r, 'lags', lags, 'peakValue', NaN, 'peakLag', NaN, 'fftpeak', NaN);
        fig = [];
        return;
    end

    % Choose the first significant peak (closest to zero lag)
    [~, iFirst] = min(lSearch(locs)); 
    bestLag     = lSearch(locs(iFirst));
    bestPk      = pks(iFirst);

    T0 = bestLag / fs;
    f0 = 1 / T0;

    % --- FFT (Welch) cross-check over the same band ---
        fchk = opts.FreqRange;
        % params for interpolation
%         df = 1; % frequency spacing
%         freq = flimits(1):df:flimits(2); % frequency interpolation axis
        [pxx,f] = welchPSD(x,fs,1,[]);
        
        iBand   = f>=fchk(1) & f<=fchk(2);
        [~,iPeak]   = max(pxx(iBand));
        f_fft   = f(iBand);
        f_fft   = f_fft(iPeak);
%         fprintf('Welch peak (cross-check): %.2f Hz\n', f_fft);

    if plotnum
        % --- Quick plots ---
        fig = figure('Color','w'); 
        subplot(2,1,1);
        lags_s = lags / fs;
        plot(lags_s, r, 'LineWidth', 1.5); hold on;
        yline(0,'k:'); 
        plot(T0, bestPk, 'ro', 'MarkerFaceColor','r');
        xlim([0, min(1/fchk(1), max(lags_s))]);
        xlabel('Lag (s)'); ylabel('Autocorr');
        title(sprintf('Autocorrelation (peak @ %.2f Hz)', f0));
    
        subplot(2,1,2);
        plot(f, 10*log10(pxx), 'LineWidth', 1.5); hold on;
        xline(f0,'r--','ACF f_0'); xline(f_fft,'k--','Welch');
        xlim(fchk); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
        title('Welch spectrum with ACF and Welch peaks');
    else 
        fig = [];
    end

    info = struct();
    info.acf        = r;
    info.lags       = lags;
    info.peakLag    = bestLag;
    info.peakValue  = bestPk;
    info.allPeaks   = table(reshape(lSearch(locs)/fs,[],1), reshape(pks,[],1), 'VariableNames', {'Lag_s','Peak'});
    info.params     = opts;
    info.fftpeak = f_fft;
end
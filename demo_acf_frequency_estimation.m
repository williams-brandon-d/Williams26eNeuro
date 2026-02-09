function demo_acf_frequency_estimation()
    % --- Simulated data: 40 Hz gamma nested in weak 8 Hz theta + noise ---
    fs  = 1000;                 % Hz
    T   = 4;                    % seconds
    t   = (0:1/fs:T-1/fs)';     
    x   = 0.7*sin(2*pi*40*t) + 0.25*sin(2*pi*8*t) + 0.6*randn(size(t));
    plotnum = 1;

    % --- Estimate with autocorrelation (restrict to, e.g., 20–90 Hz) ---
    [f0, T0, info] = acf_freq_estimate(x, fs, plotnum, ...
        'FreqRange', [20 90], ...        % Hz (set empty [] to disable)
        'Bandpass',  [20 90], ...        % Hz (set [] to skip)
        'Detrend',   true, ...
        'MaxLag',    0.5, ...            % seconds of ACF to compute
        'MinPeakProminence', 0.02);

    fprintf('ACF estimate: f0 = %.2f Hz (T0 = %.3f s), peak r = %.3f\n', f0, T0, info.peakValue);

end
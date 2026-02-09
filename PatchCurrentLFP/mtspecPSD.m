function [S,f] = mtspecPSD(data,Fs,ntapers,freq)
% params for mtftt
TW = (ntapers+1)/2; % time-bandwidth product determines number of tapers

params.Fs = Fs; % sampling frequency
params.tapers = [TW,ntapers];
params.pad = -1; % no padding
params.trialave = 0; % no trial averaging

[S,f] = mtspectrumc( data, params );

% interpolate spectra to the new frequency axis
if ~isempty(freq)
    method = 'cubic'; % interpolation method
    S = interp1(f, S, freq, method);
    f = freq;
end

end
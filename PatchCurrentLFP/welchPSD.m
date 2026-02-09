function [PSD,freq] = welchPSD(data,Fs,nSegments,freq)
% PSD parameters
% df = 1;
% f_band = [0 0.5*Fs];
% f = f_band(1):df:f_band(2);

N = length(data);

% nwin = floor(N/nSegments); % window length

ov = 0.5; % fraction window overlap
nwin = floor(N/(nSegments-(nSegments-1)*ov));

if nSegments == 1
    noverlap = 0; % no window overlap
else
    noverlap = floor(ov*nwin); % half window overlap    
end

% freq = [];

if ~isempty(freq)
    PSD = pwelch(data,nwin,noverlap,freq,Fs);
else
    pad = 0;
    nfft = max(256,2^(nextpow2(nwin)+pad));
    [PSD,freq] = pwelch(data,nwin,noverlap,nfft,Fs);
end

% figure; plot(freq,PSD); xlim([50 200]);

end
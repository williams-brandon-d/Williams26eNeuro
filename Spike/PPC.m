function ppc0 = PPC(phases)
% Vinck's Pairwise Phase Consistency (PPC) measure
% Vinck, M., Van Wingerden, M., Womelsdorf, T., Fries, P., Pennartz, C., 2010. The pairwise phase consistency: a bias-free measure of rhythmic neuronal synchronization. Neuroimage 51, 112–122.
% identical to the Aydore's unbiased PLV estimator
% Aydore, S., Pantazis, D. & Leahy, R. M., 2013. A note on the phase locking value and its properties. NeuroImage 74, 231–244.
% compute the unbiased phase locking value between spikes and a particular narrow band of the local field potential 
% the simple PLV estimate is biased for low N (# of spikes)
% the expected value of PLV > PLV for low N

% INPUT:
%  - phases: Vector of phases (radians)
% OUTPUT:
%  - ppc0: pairwise phase consistency

N = numel(phases);
scale_factor = 2/(N*(N-1));

sphases = double(sin(phases)); %% Force to double or you will rapidly get overflow
cphases = double(cos(phases));

ppc0 = 0;
for n=1:(N-1)
    ppc0 = ppc0 + ( scale_factor * ( sum( (cphases(n) .* cphases((n+1):N) ) + ( sphases(n).*sphases((n+1):N)) ) ) );
end
   
end
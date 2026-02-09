function [VS, mu, z] = vector_strength(phi)
%   VECTOR_STRENGTH Compute vector strength (mean resultant length).
%   [VS, MU] = VECTOR_STRENGTH(PHI)      % PHI in radians

    % mask NaNs
    m = ~isnan(phi);
    phi = phi(m);

    if isempty(phi)
        VS = NaN; mu = NaN; return;
    end

    % unit phasors
    z = mean(exp(1i*phi(:)));
    VS = abs(z);           % [0..1]
    mu = angle(z);         % circular mean angle in radians (-pi..pi)
end
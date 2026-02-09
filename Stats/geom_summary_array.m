function out = geom_summary_array(X, varargin)
% geom_summary_array Compute log-scale summary stats per column (default).
%
% Usage:
%   out = geom_summary_array(X)
%   out = geom_summary_array(X, 'Conf', 0.95, 'LogBase', 10, 'Dim', 1, 'VarNames', [])
%
% Inputs:
%   X         - numeric array. By default, stats are computed along Dim=1 (columns).
%
% Name-Value Pairs:
%   'Conf'    - confidence level (default 0.95)
%   'LogBase' - log base: 10 (default), 'e' (or 'ln'), or a positive numeric base ~= 1
%   'Dim'     - dimension to operate along (default 1). E.g., 1 = by column, 2 = by row
%   'VarNames'- cellstr or string array of variable names for the table output
%
% Outputs (struct fields are vectors matching the operated dimension):
%   out.n                    - sample size used per variable
%   out.geometric_mean       - geometric mean in original units
%   out.ci_lower             - lower 95% CI bound (original units)
%   out.ci_upper             - upper 95% CI bound (original units)
%   out.geometric_sd         - geometric SD (multiplicative)
%   out.mult_lo_factor       - lower multiplicative CI factor
%   out.mult_hi_factor       - upper multiplicative CI factor
%   out.has_nonpositive      - logical flag if any nonpositive values occurred in that slice
%   out.Table                - table with the above metrics
%
% Notes:
% - Requires strictly positive values in each slice (after removing NaNs) for log methods.
% - If a slice has n < 2, CI-related outputs are NaN for that slice.

% --------- Parse inputs ---------
p = inputParser;
addParameter(p, 'Conf', 0.95, @(v) isnumeric(v) && isscalar(v) && v>0 && v<1);
addParameter(p, 'LogBase', 10, @(v) (ischar(v) || isstring(v)) || (isnumeric(v) && isscalar(v) && v>0 && v~=1));
addParameter(p, 'Dim', 1, @(v) isnumeric(v) && isscalar(v) && v>=1 && v==floor(v));
addParameter(p, 'VarNames', [], @(v) isstring(v) || iscellstr(v) || isempty(v));
parse(p, varargin{:});

conf    = p.Results.Conf;
logbase = p.Results.LogBase;
dim     = p.Results.Dim;
varNames = p.Results.VarNames;

% normalize base selection
if ischar(logbase) || isstring(logbase)
    logbase = lower(string(logbase));
    if logbase=="e" || logbase=="ln"
        base = exp(1);
    else
        error('Unrecognized LogBase "%s". Use 10, ''e''/''ln'', or a positive numeric base ~= 1.', logbase);
    end
else
    base = logbase; % numeric base
end

% define log/exp functions for chosen base
if abs(base - exp(1)) < 1e-12
    logfun = @(v) log(v);
    expfun = @(v) exp(v);
else
    lb = base;
    logfun = @(v) log(v) ./ log(lb); % change-of-base
    expfun = @(v) lb.^v;
end

% --------- Prep along dimension ---------
% Permute so we operate along 1st dim
nd = ndims(X);
perm = 1:nd;
perm([1, dim]) = perm([dim, 1]);
Xp = permute(X, perm);
sz = size(Xp);
nSlices = prod(sz(2:end));

% reshape to 2D: [len x nSlices]
Xp = reshape(Xp, sz(1), nSlices);

% allocate outputs
n_vec       = nan(1, nSlices);
gm_vec      = nan(1, nSlices);
ci_lo_vec   = nan(1, nSlices);
ci_hi_vec   = nan(1, nSlices);
gsd_vec     = nan(1, nSlices);
mlo_vec     = nan(1, nSlices);
mhi_vec     = nan(1, nSlices);
nonpos_flag = false(1, nSlices);

alpha = 1 - conf;

% --------- Loop over slices (columns) ---------
for j = 1:nSlices
    x = Xp(:, j);
    x = x(isfinite(x));        % drop NaN/Inf
    if isempty(x)
        continue;              % leave as NaN
    end
    nonpos_flag(j) = any(x <= 0);
    if nonpos_flag(j)
        % cannot compute log-based stats → leave as NaN
        continue;
    end

    lx = logfun(x);
    n  = numel(lx);
    n_vec(j) = n;

    m  = mean(lx);
    s  = std(lx, 0);          % sample SD
    se = s / sqrt(max(n,1));

    % back-transform center
    gm = expfun(m);
    gm_vec(j) = gm;

    % geometric SD (multiplicative)
    gsd_vec(j) = expfun(s);

    if n >= 2
        tcrit = tinv(1 - alpha/2, n - 1);
        ci_lo_vec(j) = expfun(m - tcrit*se);
        ci_hi_vec(j) = expfun(m + tcrit*se);
        mlo_vec(j)   = expfun(-tcrit*se);  % multiplicative factors
        mhi_vec(j)   = expfun( tcrit*se);
    else
        % insufficient df for CI
        ci_lo_vec(j) = NaN;
        ci_hi_vec(j) = NaN;
        mlo_vec(j)   = NaN;
        mhi_vec(j)   = NaN;
    end
end

% --------- Pack outputs ---------
out = struct();
out.n               = n_vec;
out.geometric_mean  = gm_vec;
out.ci_lower        = ci_lo_vec;
out.ci_upper        = ci_hi_vec;
out.geometric_sd    = gsd_vec;
out.mult_lo_factor  = mlo_vec;
out.mult_hi_factor  = mhi_vec;
out.has_nonpositive = nonpos_flag;

% shape outputs back to operated dimension
outSize = sz(2:end);
fn = fieldnames(out);
for k = 1:numel(fn)
    v = out.(fn{k});
    out.(fn{k}) = reshape(v, outSize,[]);
end

% Build table for 2D use-cases (common scenario)
if nd == 2
    pVars = size(X, setdiff(1:2, dim)); % number of variables (slices)
    % default var names
    if isempty(varNames)
        varNames = "Var" + string(1:pVars);
    else
        varNames = string(varNames);
        if numel(varNames) ~= pVars
            error('VarNames length (%d) must match the number of slices (%d).', numel(varNames), pVars);
        end
    end

    % orient as row vector to match columns for Dim=1 (by column)
    tbl = table( ...
        out.n(:), ...
        out.geometric_mean(:), ...
        out.ci_lower(:), ...
        out.ci_upper(:), ...
        out.geometric_sd(:), ...
        out.mult_lo_factor(:), ...
        out.mult_hi_factor(:), ...
        out.has_nonpositive(:), ...
        'VariableNames', {'n','GM','CI_Lower','CI_Upper','GSD','Mult_CI_LowerFactor','Mult_CI_UpperFactor','HasNonPositive'} ...
    );

    % If operating along columns (Dim=1), each column is a variable (row per var)
    % If operating along rows (Dim=2), each row is a variable (row per var)
    tbl.Properties.RowNames = cellstr(varNames(:));
    out.Table = tbl;
else
    out.Table = []; % not constructed for ND arrays
end

end
